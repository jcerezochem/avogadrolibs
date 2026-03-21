/****************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the 3-Clause BSD License, (see "LICENSE").
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/core/constraint.h>

using Avogadro::Core::Constraint;
using Avogadro::MaxIndex;

TEST(ConstraintTest, ScanMetadataDefaultsToNoScan)
{
  Constraint constraint(0, 1, MaxIndex, MaxIndex, 1.5);

  EXPECT_FALSE(constraint.hasScan());
  EXPECT_EQ(constraint.scanSteps(), 0);
  EXPECT_DOUBLE_EQ(constraint.scanInitial(), 0.0);
  EXPECT_DOUBLE_EQ(constraint.scanEnd(), 0.0);
}

TEST(ConstraintTest, SetAndClearScanMetadata)
{
  Constraint constraint(0, 1, MaxIndex, MaxIndex, 1.5);
  constraint.setScan(1.2, 1.8, 7);

  EXPECT_TRUE(constraint.hasScan());
  EXPECT_DOUBLE_EQ(constraint.scanInitial(), 1.2);
  EXPECT_DOUBLE_EQ(constraint.scanEnd(), 1.8);
  EXPECT_EQ(constraint.scanSteps(), 7);

  constraint.clearScan();

  EXPECT_FALSE(constraint.hasScan());
  EXPECT_EQ(constraint.scanSteps(), 0);
  EXPECT_DOUBLE_EQ(constraint.scanInitial(), 0.0);
  EXPECT_DOUBLE_EQ(constraint.scanEnd(), 0.0);
}

TEST(ConstraintTest, CopyAndResetScanMetadata)
{
  Constraint constraint(0, 1, 2, 3, 180.0);
  constraint.setScan(-180.0, 180.0, 13);

  Constraint copied = constraint;
  EXPECT_TRUE(copied.hasScan());
  EXPECT_DOUBLE_EQ(copied.scanInitial(), -180.0);
  EXPECT_DOUBLE_EQ(copied.scanEnd(), 180.0);
  EXPECT_EQ(copied.scanSteps(), 13);

  copied.set(0, 1, 2, 3, 60.0);
  EXPECT_FALSE(copied.hasScan());
  EXPECT_EQ(copied.scanSteps(), 0);
  EXPECT_DOUBLE_EQ(copied.value(), 60.0);
}
