#include <gtest/gtest.h>

// #include "konstantinov_s_graham/all/include/ops_all.hpp"
#include "konstantinov_s_graham/common/include/common.hpp"
// #include "konstantinov_s_graham/omp/include/ops_omp.hpp"
#include "konstantinov_s_graham/seq/include/ops_seq.hpp"
// #include "konstantinov_s_graham/stl/include/ops_stl.hpp"
// #include "konstantinov_s_graham/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace konstantinov_s_graham {

class KonstantinovSRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType test_input_;
  OutType test_expected_output_;

  void SetUp() override {
    const unsigned long long duplcount = 10000000;
    const unsigned long long seedcount = 5;
    double seedx[] = {-0.1, -0.2, -0.05, 0.1, 0.05};
    double seedy[] = {-0.1, 0.1, 0.3, 0.1, -0.3};  // -0.1_-0.1 -0.2_0.1 -0.05_0.3 0.1_0.1 0.05_-0.3
    double mult = 1.0;
    test_input_.first.resize(seedcount * duplcount);
    test_input_.second.resize(seedcount * duplcount);
    for (int i = 0; i < duplcount; i++) {
      for (int j = 0; j < seedcount; j++) {
        test_input_.first[j] = seedx[j] * mult;
        test_input_.second[j] = seedy[j] * mult;
      }
      mult += 1;
    }
    test_expected_output_.resize(seedcount);
    for (int j = 0; j < seedcount; j++) {
      test_expected_output_[j] = {seedx[j] * mult, seedy[j] * mult};
      // std::cout<<seedx[j]*mult << " "<< seedy[j]*mult<<"\n";
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == test_expected_output_;
  }

  InType GetTestInputData() final {
    return test_input_;
  }
};

TEST_P(KonstantinovSRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KonstantinovAGrahamSEQ>(PPC_SETTINGS_konstantinov_s_graham);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KonstantinovSRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KonstantinovSRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace konstantinov_s_graham
