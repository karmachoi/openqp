#ifdef OQP_LIBINT_CXX_ENGINE

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <mutex>
#include <vector>

#include <libint2.hpp>

namespace {

void ensure_libint_initialized() {
  static std::once_flag init_flag;
  std::call_once(init_flag, []() { libint2::initialize(); });
}

libint2::Shell make_shell(const int32_t max_nprim, const int32_t am,
                          const int32_t pure, const int32_t nprim,
                          const double* exps, const double* coeffs,
                          const double* center) {
  libint2::svector<double> alpha;
  libint2::svector<double> c;
  alpha.reserve(nprim);
  c.reserve(nprim);
  for (int32_t p = 0; p < nprim; ++p) {
    alpha.push_back(exps[p]);
    c.push_back(coeffs[p]);
  }

  return libint2::Shell{std::move(alpha),
                        {{am, pure != 0, std::move(c)}},
                        {{center[0], center[1], center[2]}},
                        false};
}

}  // namespace

extern "C" int oqp_libint_engine_eri0(
    int32_t max_nprim, const int32_t* am, const int32_t* pure,
    const int32_t* nprim, const double* exps, const double* coeffs,
    const double* centers, double* out, int32_t out_size) {
  try {
    ensure_libint_initialized();

    std::array<libint2::Shell, 4> shells;
    auto engine_max_nprim = int32_t{0};
    auto max_l = int32_t{0};
    auto expected_size = int32_t{1};
    for (int32_t s = 0; s < 4; ++s) {
      engine_max_nprim = std::max(engine_max_nprim, nprim[s]);
      max_l = std::max(max_l, am[s]);
      shells[s] = make_shell(max_nprim, am[s], pure[s], nprim[s],
                             exps + s * max_nprim,
                             coeffs + s * max_nprim,
                             centers + s * 3);
      expected_size *= static_cast<int32_t>(shells[s].size());
    }

    if (out_size < expected_size) return 2;

    libint2::Engine engine(libint2::Operator::coulomb, engine_max_nprim, max_l, 0);
    engine.compute(shells[0], shells[1], shells[2], shells[3]);
    const auto* buf = engine.results()[0];
    if (buf == nullptr) {
      std::fill(out, out + expected_size, 0.0);
      return 1;
    }

    std::copy(buf, buf + expected_size, out);
    return 0;
  } catch (const std::exception&) {
    return -1;
  } catch (...) {
    return -2;
  }
}

#endif  // OQP_LIBINT_CXX_ENGINE
