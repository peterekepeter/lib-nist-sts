#pragma once

/* --------------------------------------------------------------------------

The following code is distributed under the following BSD-style license:

Copyright � 2013-2014 Marek Sys (syso@fi.muni.cz) & Zdenek Riha (zriha@fi.muni.cz).
All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

3. The name of the author may not be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AUTHORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------- */

inline size_t get_mask(size_t size) {
	return (1 << size) - 1;
}

inline unsigned int get_nth_block4(const unsigned char* arr, const int offset)
{
	return (*reinterpret_cast<const unsigned int*>(arr + (offset >> 3))) >> (offset & 7);//(array2[(offset >> 3)&3][(offset >> 3)] >> (offset & 7));
}

inline unsigned int get_nth_block_effect(const unsigned char* arr, const int offset)
{
	int shift = (offset & 7);
	int byte = (offset >> 3);
	if (shift == 0) return (*(unsigned int*)(arr + byte) >> shift);
	else return (*reinterpret_cast<const unsigned int*>(arr + byte) >> shift) ^ (*(unsigned int*)(arr + byte + 4) << (32 - shift));
}

inline int Mirrored_int(unsigned int val, int m) {
	int res = 0, i;
	for (i = 0; i < m; i++)
	{
		if (val & (1 << i)) res += (1 << (m - 1 - i));
	}
	return res;
}

inline void Histogram(int bitstart, int* P, int m, int bitend, const unsigned char* array) {
	int help, mask, i;
	const unsigned char* pbyte;

	mask = (1 << m) - 1;

	i = bitstart;
	while (i % 8 != 0 && i <  bitend - m + 1) {
		help = get_nth_block4(array, i);
		++P[help & mask];
		i++;
	}

	pbyte = array + i / 8;
	help = get_nth_block4(array, i);

	for (; i < bitend - m + 1 - 8; i += 8) {
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;

		++P[help & mask];
		help = *(unsigned int*)(++pbyte);
	}

	for (; i < bitend - m + 1; i++) {
		help = get_nth_block4(array, i);
		++P[help & mask];
	}
}
