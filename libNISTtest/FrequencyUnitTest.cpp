#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../libNIST/BitSequence.h"
#include "../libNIST/Parameters.h"
#include "../libNIST/Test.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace libNISTtest
{
	TEST_CLASS(FrequencyUnitTest)
	{
	public:

		TEST_METHOD(FrequencySlow)
		{
			std::vector<bool> data(8);
			
			data[0] = true;
			data[1] = false;
			data[2] = true;
			data[3] = false;
			data[4] = true;
			data[5] = false;
			data[6] = true;
			data[7] = false;

			Nist::BitSequence sequence(data);
			Nist::Parameters paramsA;

			paramsA.fast = 0;
			Nist::Test test(&sequence, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");

			// TODO: Your test code here
		}


		TEST_METHOD(FrequencyFast)
		{
			std::vector<bool> data(8);

			data[0] = true;
			data[1] = false;
			data[2] = true;
			data[3] = false;
			data[4] = true;
			data[5] = false;
			data[6] = true;
			data[7] = false;

			Nist::BitSequence sequence(data);
			Nist::Parameters paramsA;

			paramsA.fast = 1;
			Nist::Test test(&sequence, &paramsA);

			auto results = test.RunFrequency();
			Assert::IsTrue(results.sum == 0.0, L"Sum should be 0");

			// TODO: Your test code here
		}

	};
}