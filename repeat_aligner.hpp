#pragma once

#include "chaining.hpp"

#include <functional>
#include <string>
#include <string_view>

std::string repeat_aligner(const std::string &query, const std::string &target);

void call_mums(std::string_view query,
               std::string_view query_rc,
               std::string_view target,
               size_t min_length,
               const std::function<void(Ranges&&)> &callback,
               SOffset qbegin = 0,
               SOffset qend = -1,
               SOffset qrcbegin = 0,
               SOffset qrcend = -1,
               SOffset tbegin = 0,
               SOffset tend = -1);