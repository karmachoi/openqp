#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <libint2.hpp>
#include <libint2/solidharmonics.h>

namespace {

libint2::Shell make_d_shell(bool pure) {
  return libint2::Shell{{1.2, 0.4},
                        {{2, pure, {0.7, 0.3}}},
                        {{0.1, -0.2, 0.3}}};
}

std::vector<double> transform_dddd_to_pure(const double* cart) {
  constexpr auto l = 2;
  constexpr auto ncart = 6;
  constexpr auto npure = 5;
  std::vector<double> scratch_a(npure * ncart * ncart * ncart);
  std::vector<double> scratch_b(npure * npure * ncart * ncart);
  std::vector<double> scratch_c(npure * npure * npure * ncart);
  std::vector<double> pure(npure * npure * npure * npure);

  libint2::solidharmonics::transform_first(l, ncart * ncart * ncart, cart,
                                           scratch_a.data());
  libint2::solidharmonics::transform_inner(npure, l, ncart * ncart,
                                           scratch_a.data(), scratch_b.data());
  libint2::solidharmonics::transform_inner(npure * npure, l, ncart,
                                           scratch_b.data(), scratch_c.data());
  libint2::solidharmonics::transform_last(npure * npure * npure, l,
                                          scratch_c.data(), pure.data());
  return pure;
}

}  // namespace

int main() {
  libint2::initialize();

  const auto cart_shell = make_d_shell(false);
  const auto pure_shell = make_d_shell(true);

  libint2::Engine engine(libint2::Operator::coulomb, cart_shell.nprim(), 2, 0);

  engine.compute(cart_shell, cart_shell, cart_shell, cart_shell);
  const auto* cart = engine.results()[0];
  if (cart == nullptr) {
    std::cerr << "Cartesian ERI buffer is null\n";
    return 2;
  }
  std::vector<double> cart_copy(cart, cart + 6 * 6 * 6 * 6);

  engine.compute(pure_shell, pure_shell, pure_shell, pure_shell);
  const auto* pure = engine.results()[0];
  if (pure == nullptr) {
    std::cerr << "Pure ERI buffer is null\n";
    return 3;
  }
  std::vector<double> pure_copy(pure, pure + 5 * 5 * 5 * 5);

  const auto transformed = transform_dddd_to_pure(cart_copy.data());
  double max_abs_diff = 0.0;
  for (std::size_t i = 0; i < transformed.size(); ++i) {
    max_abs_diff = std::max(max_abs_diff, std::abs(transformed[i] - pure_copy[i]));
  }

  std::cout << "cartesian_size=" << cart_shell.size() << "^4\n";
  std::cout << "pure_size=" << pure_shell.size() << "^4\n";
  std::cout << "max_abs_diff=" << max_abs_diff << "\n";

  libint2::finalize();
  return max_abs_diff < 1.0e-10 ? 0 : 4;
}
