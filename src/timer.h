#pragma once
#include <iostream>
#include <time.h>

#include <iostream>
#include <string>
#include <chrono>
 
class Timer {
public:
    Timer() : _name("Time elapsed:")
    {
        restart();
    }

    explicit Timer(const std::string &name) : _name(name + ":")
    {
        restart();
    }

    inline void Restart()
    {
        restart();
    }

    inline void Restart(const std::string &name)
    {
        _name = name + ":";
        restart();
    }

    void Report(const std::string &tip, bool is_ms = true)
    {
        if (is_ms)
        {
            if (tip.length() > 0)
                std::cout << tip + ":" << elapsed() << "ms" << std::endl;
            else
                std::cout << _name << elapsed() << "ms" << std::endl;
        }
        else
        {
            if (tip.length() > 0)
                std::cout << tip + ":" << elapsed() / 1000.0 << "s" << std::endl;
            else
                std::cout << _name << elapsed() / 1000.0 << "s" << std::endl;
        }
    }
 
private:
    inline void restart()
    {
        start_time_ = std::chrono::steady_clock::now();
    }

    inline double elapsed()
    {
        end_time_ = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = end_time_ - start_time_;
        return diff.count() * 1000;
    }

private:
    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point end_time_;
    std::string _name;
}; // 