#pragma once

#include <deque>

#include <cassert>
#include <cstddef>

template <typename T, typename Container = std::deque<T>>
class OffsetVector {
  public:
    using const_iterator = typename Container::const_iterator;

    template <typename... Args>
    OffsetVector(Args&&... args) : data_(std::forward<Args>(args)...), offset_(0) {}

    void clear() {
        data_.clear();
        offset_;
    }

    template <typename... Args>
    T& emplace_back(Args&&... args) {
        return data_.emplace_back(std::forward<Args>(args)...);
    }

    T& operator[](size_t i) {
        assert(i < size());
        assert(i >= offset_);
        assert(i - offset_ < data_.size());
        return data_[i - offset_];
    }

    const T& operator[](size_t i) const {
        assert(i < size());
        assert(i >= offset_);
        assert(i - offset_ < data_.size());
        return data_[i - offset_];
    }

    T& back() { return data_.back(); }
    const T& back() const { return data_.back(); }

    T& front() { return data_.front(); }
    const T& front() const { return data_.front(); }

    const_iterator begin() const { return data_.begin(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cend() const { return data_.cend(); }

    void dealloc_front(size_t num = 1) {
        assert(data_.size() >= num);
        if (num == 1) {
            data_.pop_front();
            ++offset_;
        } else {
            data_.erase(data_.begin(), data_.begin() + num);
            offset_ += num;
        }
    }

    size_t size() const { return data_.size() + offset_; }
    size_t alloc_size() const { return data_.size(); }
    size_t offset() const { return offset_; }

    void resize(size_t new_size) {
        if (new_size < size())
            throw std::runtime_error("Shrinking not implemented");

        assert(new_size > offset_);
        data_.resize(new_size - offset_);
        assert(size() == new_size);
    }

    Container& data() { return data_; }
    const Container& data() const { return data_; }

  private:
    Container data_;
    size_t offset_;
};