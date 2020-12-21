#ifndef PROTONTRACKRATIOS
#define PROTONTRACKRATIOS

#include <utility>
#include <map>
#include <string>
/*
Helper used to calculate real 0 proton tag per arm using numbers extracted by J. Kaspar from 2017 data
https://indico.cern.ch/event/935869/contributions/3932180/attachments/2070099/3474989/kaspar_prob.pdf

Code developed by: Enrico Robutti
*/

using namespace std;

struct eraXangle
{
  pair<int, int> runRange;
  float xAngle;
  eraXangle(int firstRun, int lastRun, float angle) : runRange(firstRun, lastRun), xAngle(angle) {};
  bool operator==(const eraXangle& o) const {
    return runRange.first == o.runRange.first &&
           runRange.second == o.runRange.second &&
           xAngle == o.xAngle;};
  bool operator<(const eraXangle& o) const {
    return runRange.first < o.runRange.first ||
           (runRange.first == o.runRange.first && runRange.second < o.runRange.second) ||
           (runRange.first == o.runRange.first && runRange.second == o.runRange.second && xAngle < o.xAngle);};
  bool contains (int run, float angle) const {return (run >= runRange.first && run <= runRange.second && angle == xAngle);}
};

class protonTrackRatios
{
 public:
  protonTrackRatios() {};
  protonTrackRatios(string dataFileName) {readFromFile(dataFileName);};
  ~protonTrackRatios() {};
  int readFromFile(string dataFileName);
  float twoTracksRatio(int run, float xAngle);
  float zeroTracksRatio(int run, float xAngle, int arm);
  float trueZeroTracksRatio(int run, float xAngle, int arm);
  bool isTrueZeroTrack(int run, float xAngle, int arm);

 private:
  map<eraXangle, float > mTwoTracksRatio_;
  map<eraXangle, pair<float, float> > mZeroTracksRatio_;
  map<eraXangle, pair<float, float> > mTrueZeroTracksRatio_;
};

#endif  // PROTONTRACKRATIOS
