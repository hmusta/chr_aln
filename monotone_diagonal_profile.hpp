#pragma once

#include "helpers.hpp"

#include <algorithm>
#include <optional>
#include <vector>

#include <cassert>
#include <cstddef>

namespace chr_aln {

template <typename T>
struct Record {
    SOffset q1;
    SOffset q2;
    T penalty;

    template <typename... Args>
    Record(SOffset qq1, SOffset qq2, Args&&... ppenalty)
        : q1(qq1), q2(qq2), penalty(std::forward<Args>(ppenalty)...) {}
};

template <typename T = Penalty, class Less = std::less<T>, bool Ascending = true>
class MonotoneDiagonalProfile {
  public:
    using value_type = T;
    using record_type = Record<T>;
    using Storage = std::vector<record_type>;
    using iterator = typename Storage::iterator;
    using const_iterator = typename Storage::const_iterator;

    void clear() { records_.clear(); }
    bool empty() const { return records_.empty(); }
    size_t size() const { return records_.size(); }

    iterator begin() { return records_.begin(); }
    iterator end() { return records_.end(); }
    const_iterator begin() const { return records_.begin(); }
    const_iterator end() const { return records_.end(); }
    const_iterator cbegin() const { return records_.cbegin(); }
    const_iterator cend() const { return records_.cend(); }

    template <typename... Args>
    bool extend(SOffset q1, SOffset q2, Args&&... penalty_args) {
        if (q1 >= q2)
            return false;

        value_type penalty(std::forward<Args>(penalty_args)...);

        if (!records_.empty()) {
            iterator last = end() - 1;
            assert(!Less()(penalty, last->penalty) && "penalties must be fed non-decreasing");

            if constexpr (Ascending) {
                assert((records_.size() < 2 || records_[records_.size() - 2].q2 <= last->q1) && "ensure that stored Records are disjoint");
            } else {
                assert((records_.size() < 2 || last->q2 <= records_[records_.size() - 2].q1) && "ensure that stored Records are disjoint");
            }

            if (penalty == last->penalty) {
                bool extended = false;
                if constexpr (Ascending) {
                    if (q2 > last->q2) {
                        last->q2 = q2;
                        extended = true;
                    }
                } else {
                    if (q1 < last->q1) {
                        last->q1 = q1;
                        extended = true;
                    }
                }

                return extended;
            } else {
                if constexpr (Ascending) {
                    q1 = std::max(q1, last->q2);
                } else {
                    q2 = std::min(q2, last->q1);
                }

                if (q1 >= q2)
                    return false;
            }
        }

        if constexpr (Ascending) {
            assert((records_.empty() || records_.back().q2 <= q1) && "ensure that appended Records are disjoint");
        } else {
            assert((records_.empty() || q2 <= records_.back().q1) && "ensure that appended Records are disjoint");
        }

        records_.emplace_back(q1, q2, std::move(penalty));
        return true;
    }

    std::pair<const_iterator, const_iterator> query_min_range(SOffset a, SOffset b) const {
        if (a >= b || records_.empty())
            return std::make_pair(records_.cend(), records_.cend());

        Probe p = probe<true>(a, b);
        return std::make_pair(p.begin, p.end);
    }

    std::pair<const_iterator, const_iterator> query_min(SOffset a, SOffset b) const {
        if (a >= b || records_.empty())
            return std::make_pair(records_.cend(), records_.cend());

        Probe p = probe<false>(a, b);
        return std::make_pair(p.begin, p.end);
    }

    std::optional<Penalty> get_global_min() const {
        if (records_.empty())
            return std::nullopt;

        return records_.front().penalty;
    }

    SOffset get_global_begin() const {
        return Ascending ? records_.front().q1 : records_.back().q1;
    }

    SOffset get_global_end() const {
        return Ascending ? records_.back().q2 : records_.front().q2;
    }

  private:
    struct Probe {
        const_iterator begin;
        const_iterator end;
    };

    template <bool get_end = true>
    Probe probe(SOffset a, SOffset b) const {
        assert(a < b && "Out-of-bounds access to records in Probe::probe");

        if (records_.empty() || a >= get_global_end())
            return { records_.cend(), records_.cend() };

        if constexpr (Ascending) {
            // find the first record where (r.q2 <= a) is false <=> (r.q2 > a) is true
            const auto it = std::upper_bound(
                records_.cbegin(), records_.cend(), a,
                [](SOffset value, const record_type& r) { return r.q2 > value; });

            assert(it != records_.cend());

            // find the first record where (b <= r.q2) is true <=> (b > r.q2) is false
            return {
                it,
                get_end ? std::lower_bound(
                            it + 1, records_.cend(), b,
                            [](const record_type& r, SOffset value) { return value > r.q2; })
                        : records_.cend()
            };

        } else {
            // this side is more difficult since r.q1 is uncertain, meaning that
            // there may be an it s.t. it->q1 > (it + 1)->q2

            // find the first record where (b >= r.q2) is true <=> (b < r.q2) is false
            const auto it = std::lower_bound(
                records_.cbegin(), records_.cend(), b,
                [](const record_type& r, SOffset value) { return value < r.q2; });

            assert(it != records_.cend());

            if (b > it->q2 && it != records_.cbegin() && (it - 1)->q1 > it->q2)
                --it;

            // find the first record where (a >= r.q2) is true <=> (a < r.q2) is false
            return {
                it,
                get_end ? std::lower_bound(
                            it + 1, records_.cend(), a,
                            [](const record_type& r, SOffset value) { return value < r.q2; })
                        : records_.cend()
            };
        }
    }

    Storage records_;
};

using ForwardDiagonalProfile = MonotoneDiagonalProfile<Penalty, std::less<Penalty>, true>;
using BackwardDiagonalProfile = MonotoneDiagonalProfile<Penalty, std::less<Penalty>, false>;

} // namespace chr_aln