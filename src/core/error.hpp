#pragma once

#include <stdexcept>
#include <string>

namespace spodvr {

// Invalid input: names the offending parameter (input key) and explains the problem.
class InputError : public std::runtime_error {
 public:
  InputError(std::string parameter, const std::string& message)
      : std::runtime_error(parameter.empty() ? message : parameter + ": " + message),
        parameter_(std::move(parameter)) {}

  [[nodiscard]] const std::string& parameter() const { return parameter_; }

 private:
  std::string parameter_;
};

// A feature that exists in the Fortran code but is not (yet) available in this program.
class NotSupportedError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

}  // namespace spodvr
