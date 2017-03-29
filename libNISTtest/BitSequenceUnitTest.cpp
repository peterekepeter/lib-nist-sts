#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../libNIST/BitSequence.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace libNISTtest
{		
	TEST_CLASS(BitSequenceUnitTest)
	{
	public:
		
		TEST_METHOD(BitSequenceFromBoolVectorTest)
		{
			std::vector<bool> data(4);
			data[0] = true;
			data[1] = false;
			data[2] = true;
			data[3] = false;
			Nist::BitSequence sequence(data);
			Assert::IsTrue(sequence.GetN() == 4);
			Assert::IsTrue(sequence.ReadAt(0));
			Assert::IsFalse(sequence.ReadAt(1));
			Assert::IsTrue(sequence.ReadAt(2));
			Assert::IsFalse(sequence.ReadAt(3));
			// TODO: Your test code here
		}

	};
}