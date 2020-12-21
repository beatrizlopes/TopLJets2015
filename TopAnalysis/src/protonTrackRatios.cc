#include "TopLJets2015/TopAnalysis/interface/protonTrackRatios.h"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>

int protonTrackRatios::readFromFile(string dataFileName)
{
  ifstream dataFile(dataFileName);
  if (!dataFile.is_open()) {
    cout << "Error opening file " << dataFileName << "!" << endl;
    return -1;
  }
  // Skip descriptors line
  char line[256];
  dataFile.getline(line, 256);

  int firstRun, lastRun;
  float xAngle;
  float twoTracksRatio, zeroTracksRatio45, zeroTracksRatio56, trueZeroTracksRatio45, trueZeroTracksRatio56;
  int nCond = 0;
  while(dataFile >> firstRun >> lastRun >> xAngle >> twoTracksRatio >> zeroTracksRatio45 >> zeroTracksRatio56 >> trueZeroTracksRatio45 >> trueZeroTracksRatio56) {
    eraXangle cond(firstRun, lastRun, xAngle);
    mTwoTracksRatio_[cond] = twoTracksRatio;
    mZeroTracksRatio_[cond] = make_pair(zeroTracksRatio45, zeroTracksRatio56);
    mTrueZeroTracksRatio_[cond] = make_pair(trueZeroTracksRatio45, trueZeroTracksRatio56);
    nCond++;
  }

  return nCond;
}


float protonTrackRatios::twoTracksRatio(int run, float xAngle)
{
  for (auto const& twoTracksRatio : mTwoTracksRatio_) {
    if (twoTracksRatio.first.contains(run, xAngle)) {
      return twoTracksRatio.second;
    }
  }
  return 1000.;
}


float protonTrackRatios::zeroTracksRatio(int run, float xAngle, int arm)
{
  for (auto const& zeroTracksRatio : mZeroTracksRatio_) {
    if (zeroTracksRatio.first.contains(run, xAngle)) {
      if (arm == 0)
	return zeroTracksRatio.second.first;
      else if (arm == 1)
	return zeroTracksRatio.second.second;
      else return 1000.;
    }
  }
  return 1000.;
}


float protonTrackRatios::trueZeroTracksRatio(int run, float xAngle, int arm)
{
  for (auto const& trueZeroTracksRatio : mTrueZeroTracksRatio_) {
    if (trueZeroTracksRatio.first.contains(run, xAngle)) {
      if (arm == 0)
	return trueZeroTracksRatio.second.first;
      else if (arm == 1)
	return trueZeroTracksRatio.second.second;
      else return 1000.;
    }
  }
  return 1000.;
}


bool protonTrackRatios::isTrueZeroTrack(int run, float xAngle, int arm)
{
  srand(time(NULL));
  float rnd = float(rand())/RAND_MAX;
  if (rnd < trueZeroTracksRatio(run, xAngle, arm))
    return true;
  else
    return false;
}
