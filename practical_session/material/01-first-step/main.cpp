// Copyright 2025 the samurai team
// SPDX-License-Identifier:  BSD-3-Clause

#include <iostream>
#include <samurai/samurai.hpp>

int main(int argc, char* argv[])
{
    samurai::initialize("samurai setup test", argc, argv);
    SAMURAI_PARSE(argc, argv);

    samurai::finalize();
    return 0;
}
