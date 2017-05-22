#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../LibNistSts/BitSequence.h"

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

		TEST_METHOD(BitSequenceMove)
		{
			std::vector<bool> data(4);
			data[0] = true;
			data[1] = false;
			data[2] = true;
			data[3] = false;
			Nist::BitSequence sequence(data);
			Nist::BitSequence sequence2 = std::move(sequence);
			Assert::IsTrue(sequence.GetEpsilon() == nullptr, L"Original sequence should not contain any data.");
			Assert::IsTrue(sequence2.GetEpsilon() != nullptr, L"New sequence should have data.");
			Assert::IsTrue(sequence2.GetN() == 4, L"GetN() for new sequence should return 4.");
		}
	};
}