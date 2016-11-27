#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../proj_engine.h"
#include <string>
using namespace std;




TEST_CASE( "divited ", "[bign_init1]" ) {
    REQUIRE( bign_init1("1") == "1" );
    REQUIRE( bign_init1("10") == "10" );
    REQUIRE( bign_init1("453678432567895432567896543215678") == "453678432567895432567896543215678" );
}

