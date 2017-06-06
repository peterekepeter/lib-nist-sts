#include "BitSequence.h"

void Nist::BitSequence::convert_epsilon_to_array()
{
	size_t byte_size, i, bytes, j;

	bytes = n >> 3;

	// align to 16 bytes
	byte_size = ((bytes + 15) >> 4) << 4;
	if (array == nullptr)
	{
		array = new unsigned char[byte_size];
	}

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

void Nist::BitSequence::convert_array_to_epsilon()
{
	if (epsilon == nullptr)
	{
		epsilon = new unsigned char[this->GetN()];
	}
	size_t bytes = n >> 3;
	// each byte
	for (size_t i=0; i<bytes; i++)
	{
		auto byte = array[i];
		for (size_t j=0; j<8; j++)
		{
			epsilon[(i << 3) + j] = (array[i] >> j) & 1;
		}
	}
	// last remaning bits
	for (size_t j = bytes << 3; j < n; j++)
	{
		epsilon[j] = (array[bytes] >> (j - (bytes << 3))) & 1;
	}
}

void Nist::BitSequence::AllocateArrayIfNull()
{
	if (array == nullptr)
	{
		size_t bytes = n >> 3;
		size_t byte_size = ((bytes + 15) >> 4) << 4;
		array = new unsigned char[byte_size];
	}
}
