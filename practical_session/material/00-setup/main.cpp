#include <iostream>
#include <samurai/samurai.hpp>
#include <samurai/box.hpp>

int main(int argc, char* argv[])
{
    samurai::initialize("samurai setup test", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::Box<double, 1> box({0.0}, {1.0});
    std::cout << box << std::endl;
    samurai::finalize();
    return 0;
}
