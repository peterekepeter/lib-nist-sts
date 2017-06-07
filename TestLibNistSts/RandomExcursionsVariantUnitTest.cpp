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
	TEST_CLASS(RandomExcursionVariantUnitTest)
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
				bit = (lcg & 8) == 8;
				lcg = lcg * 7612341 + 41235467;
			}
			return Nist::BitSequence(data);
		}

		// compare each version
		TEST_METHOD(RandomExcursionVariantVersions)
		{
			auto seq = GetSeqLcg(4096);
			seq[0] = false;
			seq[1] = false;
			seq[2] = false;
			Nist::Parameters parameters;
			Nist::Results results1, results2;
			Nist::Test test1(&seq, &parameters, &results1);
			Nist::Test test2(&seq, &parameters, &results2);
			RandomExcursionsVariant(test1);
			RandomExcursionsVariant2(test2);
			auto& result1 = test1.GetResults().random_excursion_variant;
			auto& result2 = test2.GetResults().random_excursion_variant;
			Assert::AreEqual(result1.valid, result2.valid);
			if (result1.valid && result2.valid)
			{
				for (size_t i = 0; i < result1.p_value.size(); i++) {
					Assert::AreEqual(result1.p_value[i], result2.p_value[i], 0.01);
					Assert::AreEqual(result1.x[i], result2.x[i]);
					Assert::AreEqual(result1.count[i], result2.count[i]);
				}
			}
		}

	};
}