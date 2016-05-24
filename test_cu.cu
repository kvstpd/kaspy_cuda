
#include <stdio.h>

int vectorr(void);

extern "C" int TEST_ME(int * test)
{
    printf("test %d\n", *test);
    
    vectorr();
    
    return 0;

}