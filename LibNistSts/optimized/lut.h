#pragma once

// lookup tables

extern short int LU_byte_weight[256];
extern short int LU_byte_switches[256];
extern unsigned char LU_byte_inverted[256];

extern signed char LUT_HW_8[], LUT_Switches_8[];
extern signed char LUT_HW_16[], LUT_Switches_16[];
extern signed char LUT_Lrun_start_8[], LUT_Lrun_end_8[], LUT_Lrun_max_8[];
extern signed char LUT_Lrun_start_16[], LUT_Lrun_end_16[], LUT_Lrun_max_16[];
extern signed char LUT_Cusum_max_positiv_8[], LUT_Cusum_max_negativ_8[], LUT_Cusum_8[];
extern signed char LUT_Cusum_max_positiv_16[], LUT_Cusum_max_negativ_16[], LUT_Cusum_16[];