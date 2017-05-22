#pragma once
#include <vector>

namespace Nist
{
	class BitSequence
	{
	private:
		unsigned char* epsilon = nullptr;
		unsigned char* array = nullptr;
		size_t n = 0;

		// clear the memory used by this object
		void free() {
			array = nullptr;
			epsilon = nullptr;
			//TODO fix this
			if (array != nullptr)
			{
				delete array;
			}
			if (epsilon != nullptr)
			{
				delete epsilon;
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
			for (size_t i=0; i<n; i++)
			{
				epsilon[i] = booleanVector[i];
			}
		}

		~BitSequence()
		{
			free();
		}

		const unsigned char* GetEpsilon() 
		{
			return epsilon;
		}

		const unsigned char* GetArray()
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


		BitSequence(const BitSequence& other) = delete;

		BitSequence(BitSequence&& other) noexcept
			: epsilon(other.epsilon),
			array(other.array),
			n(other.n)
		{
			other.epsilon = nullptr;
			other.array = nullptr;
			other.n = 0;
		}

		BitSequence& operator=(const BitSequence& other) = delete;

		// move assignent
		BitSequence& operator=(BitSequence&& other) noexcept
		{
			if (this == &other)
				return *this;
			// free this if necessarry
			free();
			// steal data pointers from other
			epsilon = other.epsilon;
			array = other.array;
			n = other.n;
			// set other's pointers to null
			other.epsilon = nullptr;
			other.array = nullptr;
			other.n = 0;
			return *this;
		}

	};


}
