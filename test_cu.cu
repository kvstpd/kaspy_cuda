
#include <stdio.h>

int vectorr(void);


extern "C" int test_me_(int * test)
{
    printf("test %d\n", *test);

    vectorr();

    return 0;

}


extern "C" int TEST_ME(int * test)
{
    printf("test %d\n", *test);
    
    vectorr();
    
    return 0;

}