"""PyOQP

This module serves as the main entry point for the OQP (Open Quantum Package)
software. It handles command-line arguments, initializes the OQP environment,
and executes the appropriate computations or tests.
"""
import os
import sys
import time
import argparse
from signal import signal, SIGINT, SIG_DFL
import oqp
from oqp.utils.file_utils import dump_log
from oqp.utils.input_checker import check_input_values
from oqp.molecule import Molecule
from oqp.library.runfunc import (
   compute_energy, compute_grad, compute_nac, compute_soc, compute_geom,
   compute_nacme, compute_properties, compute_hess, compute_thermo
)


class Runner:
    """
    OQP main class for running calculations and tests.

    one can read input parameters from text file using input_file, e.g,
        OQP(project=project_name, log=log_file, input_file=a_text_file)

    or

    one can pass input parameters to OQP internally using input_dict, e.g.,
        OQP(project=project_name, log=log_file, input_dict=USER_CONFIG)
        USER_CONFIG follows the same format as OQP_CONFIG_SCHEMA

    """

    def __init__(self, project=None, input_file=None, log=None,
                 input_dict=None, silent=0):
        """
        Initialize the OQP Runner.

        Args:
            project (str): Project name.
            input_file (str): Path to the input file.
            log (str): Path to the log file.
            input_dict (dict): Dictionary containing input parameters.
            silent (int): Flag to run in silent mode (0 or 1).
        """
        start_time = time.time()

        signal(SIGINT, SIG_DFL)

        if os.path.exists(log):
            os.remove(log)

        # Initialize mol
        if input_dict:
            if isinstance(input_dict, dict):
                self.mol = Molecule(project, input_file, log, silent=silent).from_dict(input_dict)
        else:
            self.mol = Molecule(project, input_file, log, silent=silent).from_file(input_file)

        # Check input values
        check_input_values(self.mol.config)

        # Reload mol if possible
        self.reload()

        # Attach the starting time to mol
        self.mol.start_time = start_time

        dump_log(self.mol, title='', section='start')

    def run(self, test_mod=False):
        """
        Run the OQP calculation or test.

        Args:
            test_mod (bool): Flag to run in test mode.
        """
        # Set up logfile
        self.mol.data["OQP::log_filename"] = self.mol.log

        # Set up banner
        oqp.oqp_banner(self.mol)

        # Set up basis set
        oqp.library.set_basis(self.mol)

        # Get the run type from mol configuration
        run_type = self.mol.config["input"]["runtype"]

        # Define the mapping of run types to their respective functions
        run_func = {
            'energy': compute_energy,
            'grad': compute_grad,
            'nac': compute_nac,
            'nacme': compute_nacme,
            'soc': compute_soc,
            'optimize': compute_geom,
            'meci': compute_geom,
            'mecp': compute_geom,
            'mep': compute_geom,
            'ts': compute_geom,
            'hess': compute_hess,
            'thermo': compute_thermo,
            'prop': compute_properties,
        }

        # Turn off convergence check for test
        if test_mod:
            self.mol.config['tests']['exception'] = True

        # Run calculations
        run_func[run_type](self.mol)

        dump_log(self.mol, title='', section='end')

    def results(self):
        """ Collect calculation results for internal Python communication.

        Returns:
            dict: A dictionary containing calculation results.
        """
        summary = {
            'atoms': self.mol.get_atoms(),
            'system': self.mol.get_system(),
            'energy': self.mol.energies,
            'grad': self.mol.grads,
            'dcm': self.mol.dcm,
            'nac': self.mol.nac,
            'soc': self.mol.soc,
            'data': self.mol.get_data(),
        }
        return summary

    def reload(self):
        """
        Reload calculation based on the specified guess type.
        """
        guess_type = self.mol.config['guess']['type']
        guess_file = self.mol.config['guess']['file']

        if guess_type == 'json':
            self.mol.load_data()
        elif guess_type == 'auto':
            if os.path.exists(guess_file):
                self.mol.load_data()
        else:
            pass

    def back_door(self, data):
        """
        Set back door data for the molecule.

        Args:
            data: Data to be set as back door for the molecule.
        """
        self.mol.back_door = data

    def test(self):
        """
        Run tests and check results against reference data (json).
        """
        message, diff = self.mol.check_ref()

        return message, diff


def main():
    """
    Main function to handle command-line arguments and run OQP.
    """
    parser = argparse.ArgumentParser(description='OQP Runner',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', nargs='?', help='Input file')
    parser.add_argument('--run_tests',
                        metavar='path',
                        help='run tests from a specified folder or:\n'
                             '  all    - Run all tests in examples\n'
                             '  other  - Run tests in examples/other')
    parser.add_argument('--silent', action='store_true', help='run silently')
    args = parser.parse_args()

    if args.run_tests:
        report, status = run_tests(args.run_tests)
        print(report)
        sys.exit(status)

    if not args.input:
        parser.print_help()
        sys.exit(1)

    input_file = os.path.abspath(args.input)
    if not os.path.exists(input_file):
        print(f'\n   PyOQP input file {input_file} not found\n')
        sys.exit(1)

    input_path = os.path.dirname(input_file)
    project_name = os.path.basename(input_file).replace('.inp', '')
    log = f'{input_path}/{project_name}.log'

    silent = 1 if args.silent else 0

    # Initialize OQP class
    oqp_runner = Runner(project=project_name,
                        input_file=input_file,
                        log=log,
                        silent=silent)

    # Run OQP
    oqp_runner.run()
    oqp_runner.results()


def run_tests(test_path):
    """
    Run OQP tests.

    Args:
        test_path (str): Path to the test directory or 'all' or 'other'.

    Returns:
        str: Test report.
    """
    from oqp.utils.oqp_tester import OQPTester
    tester = OQPTester(base_test_dir=None,
                       output_dir='openqp_test_tmp',
                       total_cpus=None,
                       omp_threads=4)
    return tester.run(test_path), tester.status


if __name__ == "__main__":
    main()