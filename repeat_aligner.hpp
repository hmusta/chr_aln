#pragma once

#include "chaining.hpp"

#include <functional>
#include <string>
#include <string_view>

#include <bindings/cpp/WFAligner.hpp>

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

std::pair<Score, std::string>
get_alignment(wfa::WFAligner& aligner,
                const ScoreModel &score_model,
                std::string_view query,
                std::string_view target,
                bool penalty_to_score = true,
                SOffset heuristics_length_cutoff = max_offset,
                Diag min_k = min_diag,
                Diag max_k = max_diag);