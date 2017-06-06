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
		void convert_epsilon_to_array();
		void convert_array_to_epsilon();
		void AllocateArrayIfNull();

	public:

		friend class BitRef;

		// a reference to a single bit inside bitsequence
		class BitRef
		{
		private:
			const BitSequence& sequence;
			size_t position;
#ifdef _DEBUG
			bool readonly;
#endif
		public:
			BitRef(BitSequence& sequence, size_t position)
				: sequence(sequence), position(position)
#ifdef _DEBUG
			, readonly(false) 
#endif
			{ }

			BitRef(const BitSequence& sequence, size_t position)
				: sequence(sequence), position(position)
#ifdef _DEBUG
				, readonly(false)
#endif
			{ }

			operator unsigned char() const
			{
#ifdef _DEBUG
				if (position>=0 && position<sequence.n)
				{
#endif
					if (sequence.array != nullptr)
					{
						return (sequence.array[position >> 3] >> (position & 7)) & 1;
					}
					else if (sequence.epsilon != nullptr)
					{
						return sequence.epsilon[position];
					}
					else
					{
						throw std::logic_error("Sequence lacks storage.");
					}
#ifdef _DEBUG
				}
				else
				{
					throw std::out_of_range("BitRef out of range.");
				}
#endif
			}

			inline operator bool() const
			{
				return static_cast<unsigned char>(*this) != 0;
			}

			// set bit to value
			inline BitRef& operator =(unsigned char value)
			{
#ifdef _DEBUG
				if (value > 1)
				{
					throw std::invalid_argument("value should be 0 or 1");
				}
				if (readonly)
				{
					throw std::invalid_argument("sequence is read only");
				}
				if (position>=0 && position<sequence.n)
				{
#endif
				if (sequence.array != nullptr)
				{
					auto& number = sequence.array[position >>3];
					auto n = position & 7;
					auto x = value;
					number = number & ~(1 << n) | (x << n);
				}
				if (sequence.epsilon != nullptr)
				{
					sequence.epsilon[position] = value;
				}
#ifdef _DEBUG
				}
				else
				{
					throw std::out_of_range("BitRef out of range.");
				}
#endif
				return *this;
			}

			BitRef& operator =(bool value){
				this->operator=(static_cast<unsigned char>(value == true ? 1 : 0));
				return *this;
			}

			BitRef& operator =(char  value) {
				this->operator=(static_cast<unsigned char>(value));
				return *this;
			}

			BitRef& operator =(size_t value) {
				this->operator=(static_cast<unsigned char>(value));
				return *this;
			}

			BitRef& operator =(int32_t value) {
				this->operator=(static_cast<unsigned char>(value));
				return *this;
			}

			BitRef& operator =(int64_t value) {
				this->operator=(static_cast<unsigned char>(value));
				return *this;
			}

		};


		BitSequence()
		{
			
		}

		BitSequence(size_t n) : n(n)
		{
			AllocateArrayIfNull();
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
			if (epsilon == nullptr)
			{
				if (array != nullptr)
				{
					convert_array_to_epsilon();
				}
			}
			return epsilon;
		}

		const unsigned char* GetArray()
		{
			if (array == nullptr)
			{
				if (epsilon != nullptr)
				{
					convert_epsilon_to_array();
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
			return BitRef(*this, pos);
		}

		BitRef operator[](size_t position)
		{
			return BitRef(*this, position);
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
