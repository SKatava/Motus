#include <gtest/gtest.h>

#include "Motus/Profiler.h"

#include <thread>

using namespace motus;

TEST(SectionStats, TracksCountTotalMinMaxLastAvg) {
    SectionStats s;
    EXPECT_DOUBLE_EQ(s.avg(), 0.0);

    s.add(2.0);
    s.add(8.0);
    s.add(5.0);

    EXPECT_EQ(s.count, 3);
    EXPECT_DOUBLE_EQ(s.total, 15.0);
    EXPECT_DOUBLE_EQ(s.min, 2.0);
    EXPECT_DOUBLE_EQ(s.max, 8.0);
    EXPECT_DOUBLE_EQ(s.last, 5.0);
    EXPECT_DOUBLE_EQ(s.avg(), 5.0);
}

TEST(Profiler, ReportListsSectionsInFirstSeenOrder) {
    //The profiler is a process wide singleton, so use names no other test records
    auto& p = Profiler::get();
    p.record("test/zeta", 1.0);
    p.record("test/alpha", 2.0);
    p.record("test/zeta", 3.0);

    testing::internal::CaptureStdout();
    p.report();
    const std::string out = testing::internal::GetCapturedStdout();

    const auto zeta  = out.find("test/zeta");
    const auto alpha = out.find("test/alpha");
    ASSERT_NE(zeta,  std::string::npos);
    ASSERT_NE(alpha, std::string::npos);
    EXPECT_LT(zeta, alpha);
    EXPECT_EQ(out.find("test/zeta", zeta + 1), std::string::npos) << "section listed twice";
}

TEST(ScopeTimer, RecordsElapsedTime) {
    {
        ScopeTimer t("test/sleep");
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    testing::internal::CaptureStdout();
    Profiler::get().report();
    const std::string out = testing::internal::GetCapturedStdout();

    const auto line = out.find("test/sleep");
    ASSERT_NE(line, std::string::npos);

    double avg = 0, min = 0, max = 0;
    long long count = 0;
    std::sscanf(out.c_str() + line, "test/sleep %lf %lf %lf %lld", &avg, &min, &max, &count);
    EXPECT_EQ(count, 1);
    EXPECT_GE(avg, 5.0);
}
