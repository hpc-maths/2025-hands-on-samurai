#include <iostream>
#include <samurai/samurai.hpp>

int main(int argc, char* argv[])
{
    samurai::initialize("samurai setup test", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::finalize();
    return 0;
}
