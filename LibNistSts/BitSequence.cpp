#include "BitSequence.h"

void Nist::BitSequence::convert_epsilon_to_array(int n)
{
	int byte_size, i, bytes, j;

	bytes = n >> 3;

	// align to 16 bytes
	byte_size = ((bytes + 15) >> 4) << 4;
	array = new unsigned char[byte_size];

	for (i = 0; i < bytes; i++)
	{
		array[i] = 0;
		for (j = 0; j < 8; j++)
			array[i] |= epsilon[(i << 3) + j] << j;
	}
	array[i] = 0;
	for (j = bytes << 3; j < n; j++)
		array[i] |= ((epsilon[j]) << (j - (bytes << 3)));

	// fill rest with 0
	for (i=i+1;i<byte_size;i+=1)
	{
		array[i] = 0;
	}
}
