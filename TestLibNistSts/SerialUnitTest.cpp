#include "stdafx.h"
#include "CppUnitTest.h"
#include <vector>
#include "../LibNistSts/BitSequence.h"
#include "../LibNistSts/Parameters.h"
#include "../LibNistSts/Test.h"
#include "../LibNistSts/common/stat_fncs.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace libNISTtest
{
	TEST_CLASS(SerialUnitTest)
	{
	public:

		Nist::BitSequence GetSeqLcg(size_t n)
		{
			std::vector<bool> data(n);
			// make sequence of 'random'
			int lcg = 112437;
			data.resize(256);
			for (auto& bit : data)
			{
				bit = lcg & 8 == 8;
				lcg = lcg * 7612341 + 41235467;
			}
			return Nist::BitSequence(data);
		}

		// compare each version
		TEST_METHOD(SerialCompareVersions)
		{
			auto seq = GetSeqLcg(4096);
			seq[0] = false;
			seq[1] = false;
			seq[2] = false;
			Nist::Parameters parameters;
			Nist::Results results1, results2, results3;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			Nist::Test test3(&seq, &parameters, &results3);
			Serial(test1);
			Serial2(test2);
			Serial4(test3);
			auto& result1 = test1.GetResults().serial;
			auto& result2 = test2.GetResults().serial;
			auto& result3 = test3.GetResults().serial;
			Assert::AreEqual(result1.del1, result2.del1, 0.01);
			Assert::AreEqual(result1.del1, result3.del1, 0.01);
			Assert::AreEqual(result1.del2, result2.del2, 0.01);
			Assert::AreEqual(result1.del2, result3.del2, 0.01);
			Assert::AreEqual(result1.psim0, result2.psim0, 0.01);
			Assert::AreEqual(result1.psim0, result3.psim0, 0.01);
			Assert::AreEqual(result1.psim1, result2.psim1, 0.01);
			Assert::AreEqual(result1.psim1, result3.psim1, 0.01);
			Assert::AreEqual(result1.psim2, result2.psim2, 0.01);
			Assert::AreEqual(result1.psim2, result3.psim2, 0.01);
			Assert::AreEqual(result1.p_value1, result2.p_value1, 0.01);
			Assert::AreEqual(result1.p_value1, result3.p_value1, 0.01);
			Assert::AreEqual(result1.p_value2, result2.p_value2, 0.01);
			Assert::AreEqual(result1.p_value2, result3.p_value2, 0.01);
		}

	};
}