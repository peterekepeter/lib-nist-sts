#pragma once
#include <vector>

namespace Nist
{
	class BitSequence
	{
	private:
		unsigned char* epsilon = nullptr;
		unsigned char *array = nullptr;
		size_t n = 0;

		void free() {
			if (epsilon != nullptr)
			{
				delete epsilon;
				epsilon = nullptr;
			}
			if (array != nullptr)
			{
				delete array;
				array = nullptr;
			}
		}

		// convert from epsilon representation to array represetntation
		void convert_epsilon_to_array(int n);

	public:

		BitSequence()
		{
			
		}

		BitSequence(const std::vector<bool> booleanVector)
		{
			n = booleanVector.size();
			epsilon = new unsigned char[booleanVector.size()];
			for (int i=0; i<n; i++)
			{
				epsilon[i] = booleanVector[i];
			}
		}

		~BitSequence()
		{
			free();
		}

		unsigned char* GetEpsilon()
		{
			return epsilon;
		}

		unsigned char* GetArray()
		{
			if (array == nullptr)
			{
				if (epsilon != nullptr)
				{
					convert_epsilon_to_array(n);
				}
			}
			return array;
		}

		size_t GetN() const
		{
			return n;
		}

		// for testing purposes, not optimised
		bool ReadAt(size_t pos) const
		{
			if (pos>=n)
			{
				throw std::out_of_range("Given index must be smaller than N.");
			}
			if (epsilon != nullptr)
			{
				return epsilon[pos] != 0;
			}
			if (array != nullptr)
			{
				return (array[pos >> 3] >> (pos&7)) != 0;
			}
			return false;
		}
	};


}
