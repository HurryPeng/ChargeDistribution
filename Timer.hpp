// Timer.hpp
// HurryPeng

#include <chrono>

namespace HurryPeng
{

class Timer
{
    public:
        Timer();

        void restart();
        long double read();

    private:
        decltype(std::chrono::system_clock::now()) begin;
};

Timer::Timer() :begin(std::chrono::system_clock::now()) {}

void Timer::restart()
{
    begin = std::chrono::system_clock::now();
}

long double Timer::read()
{
    return (long double)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - begin).count()) / 1000000;
}

} // namespace HurryPeng
