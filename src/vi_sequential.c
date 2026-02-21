/**
 * vi_sequential.c - Algorithm 3: Sequential VI for two-product substitution
 */

#include "common.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    Instance inst = INSTANCE_P1;
    printf("CUDA Perishable Inventory VI - Sequential\n");
    printf("Instance P1: mu_a=%d mu_b=%d Qa=%d Qb=%d N=%d\n",
           inst.mu_a, inst.mu_b, inst.Qa, inst.Qb, inst.N);
    printf("(Full VI implementation in Phase 3)\n");
    return 0;
}
