#include "NetworkUtils.h"
#include <arpa/inet.h>

bool isBigEndian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};
    return bint.c[0] == 1;
}

uint64_t ntohll(uint64_t value) {
    if (isBigEndian()) {
        return value;
    }
    return ((uint64_t)ntohl(value & 0xFFFFFFFF) << 32) | ntohl(value >> 32);
}
