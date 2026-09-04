/***********************************************/
/**
* @file gnssReceiver.cpp
*
* @brief GNSS receiver.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#include <random>
#include "base/import.h"
#include "base/string.h"
#include "parser/expressionParser.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "inputOutput/logging.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnssLambda.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"

/***********************************************/

GnssReceiver::GnssReceiver(Bool isMyRank, Bool isEarthFixed, const Platform &platform,
                           GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs,
                           Bool integerAmbiguities, Double wavelengthFactor)
  : GnssTransceiver(platform, noPatternFoundAction, useableEpochs), isMyRank_(isMyRank),
    isEarthFixed_(isEarthFixed), integerAmbiguities(integerAmbiguities), wavelengthFactor(wavelengthFactor)
{
}

/***********************************************/

void GnssReceiver::copyObservations2ContinuousMemoryBlock()
{
  try
  {
    UInt count = 0;
    for(const auto &obsEpoch : observations_)
      for(const GnssObservation *obs : obsEpoch)
        if(obs)
          count++;
    obsMem.resize(count);
    obsMem.shrink_to_fit();
    // copy
    count = 0;
    for(auto &obsEpoch : observations_)
      for(GnssObservation *&obs : obsEpoch)
        if(obs)
        {
          obsMem.at(count) = *obs;
          obsMem.at(count).shrink_to_fit();
          delete obs;
          obs = &obsMem.at(count++);
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::disable(UInt idEpoch, const std::string &reason)
{
  try
  {
    GnssTransceiver::disable(idEpoch, reason);
    if(idEpoch < observations_.size())
      observations_.at(idEpoch).clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::disable(const std::string &reason)
{
  try
  {
    GnssTransceiver::disable(reason);
    if(!reason.empty() && isMyRank())
      disableReason = reason;
    isMyRank_ = FALSE;
    obsMem.clear();
    obsMem.shrink_to_fit();
    observations_.clear();
    observations_.shrink_to_fit();
    tracks.clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssObservation *GnssReceiver::observation(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < observations_.size()) && (idTrans < observations_.at(idEpoch).size()))
    return observations_[idEpoch][idTrans];
  return nullptr;
}

/***********************************************/

void GnssReceiver::deleteObservation(UInt idTrans, UInt idEpoch)
{
  try
  {
    if(!observation(idTrans, idEpoch))
      return;
    observations_[idEpoch][idTrans] = nullptr;
    if(std::all_of(observations_[idEpoch].begin(), observations_[idEpoch].end(), [](auto obs) {return obs == nullptr;}))
      disable(idEpoch, "no valid transmitters left");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::signalCompositionDefault(const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A)
{
  try
  {
    // composed type = factor1 * type1 + factor2 * type2
    static const std::vector<std::tuple<GnssType, GnssType, GnssType, Double, Double>> composites =
      {{GnssType::C1XG,  GnssType::C1SG, GnssType::C1LG, 0.5, 0.5},
       {GnssType::C2XG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5XG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},

       {GnssType::C4XR,  GnssType::C4AR, GnssType::C4BR, 0.5, 0.5},
       {GnssType::C6XR,  GnssType::C6AR, GnssType::C6BR, 0.5, 0.5},
       {GnssType::C3XR,  GnssType::C3IR, GnssType::C3QR, 0.5, 0.5},

       {GnssType::C1XE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5XE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7XE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8XE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6XE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5},

       {GnssType::C2XC,  GnssType::C2IC, GnssType::C2QC, 0.5, 0.5},
       {GnssType::C1XC,  GnssType::C1DC, GnssType::C1PC, 0.5, 0.5},
       {GnssType::C1ZC,  GnssType::C1SC, GnssType::C1LC, 0.5, 0.5},
       {GnssType::C5XC,  GnssType::C5DC, GnssType::C5PC, 0.5, 0.5},
       {GnssType::C7XC,  GnssType::C7IC, GnssType::C7QC, 0.5, 0.5},
       {GnssType::C7ZC,  GnssType::C7DC, GnssType::C7PC, 0.5, 0.5},
       {GnssType::C8XC,  GnssType::C8DC, GnssType::C8PC, 0.5, 0.5},
       {GnssType::C6XC,  GnssType::C6IC, GnssType::C6QC, 0.5, 0.5},
       {GnssType::C6ZC,  GnssType::C6DC, GnssType::C6PC, 0.5, 0.5},

       {GnssType::C1XJ,  GnssType::C1SJ, GnssType::C1LJ, 0.5, 0.5},
       {GnssType::C2XJ,  GnssType::C2SJ, GnssType::C2LJ, 0.5, 0.5},
       {GnssType::C5XJ,  GnssType::C5IJ, GnssType::C5QJ, 0.5, 0.5},
       {GnssType::C5ZJ,  GnssType::C5DJ, GnssType::C5PJ, 0.5, 0.5},
       {GnssType::C6XJ,  GnssType::C6SJ, GnssType::C6LJ, 0.5, 0.5},
       {GnssType::C6ZJ,  GnssType::C6SJ, GnssType::C6EJ, 0.5, 0.5},

       // unknown attributes
       {GnssType::C2UG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5UG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},
       {GnssType::C1UE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5UE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7UE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8UE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6UE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5}};

    typesTrans = GnssType::replaceCompositeSignals(types);

    A = Matrix(types.size(), typesTrans.size());
    for(UInt idType=0; idType<types.size(); idType++)
      if((types.at(idType) == GnssType::PHASE) || (types.at(idType) == GnssType::RANGE)) // only phase and code signals are transmitted (what about doppler?)
      {
        GnssType type = types.at(idType);

        UInt idx;
        if(type.isInList(typesTrans, idx))
        {
          A(idType, idx) = 1.; // signal observed directly
          continue;
        }

        if(type == GnssType::C2DG)
        {
          A(idType, GnssType::index(typesTrans, GnssType::C1CG)) = +1.;
          A(idType, GnssType::index(typesTrans, GnssType::C1WG)) = -1.;
          A(idType, GnssType::index(typesTrans, GnssType::C2WG)) = +1.;
          continue;
        }

        const auto composite = std::find_if(composites.begin(), composites.end(), [&](const auto &composite) {return (std::get<0>(composite) == type);});
        if(composite != composites.end())
        {
          A(idType, GnssType::index(typesTrans, std::get<1>(*composite))) = std::get<3>(*composite);
          A(idType, GnssType::index(typesTrans, std::get<2>(*composite))) = std::get<4>(*composite);
          continue;
        }

        throw(Exception("composite signal not implemented: "+type.str()));
      } // for(idType)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::signalComposition(UInt /*idEpoch*/, const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A) const
{
  signalCompositionDefault(types, typesTrans, A);
}

/***********************************************/
/***********************************************/

GnssReceiver::ObservationEquationList::ObservationEquationList(const GnssReceiver &receiver, const std::vector<GnssTransmitterPtr> &transmitters,
                                                               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                                               const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                               GnssObservation::Group group)
{
  try
  {
    eqn.resize(receiver.idEpochSize());
    for(UInt idEpoch=0; idEpoch<eqn.size(); idEpoch++)
    {
      eqn.at(idEpoch).resize(receiver.idTransmitterSize(idEpoch));
      for(UInt idTrans=0; idTrans<eqn.at(idEpoch).size(); idTrans++)
      {
        GnssObservation *obs = receiver.observation(idTrans, idEpoch);
        std::vector<GnssType> types;
        if(obs && obs->observationList(group, types))
        {
          auto e = new GnssObservationEquation(*obs, receiver, *transmitters.at(idTrans), rotationCrf2Trf, reduceModels, idEpoch, FALSE, types);
          eqn.at(idEpoch).at(idTrans) = std::unique_ptr<GnssObservationEquation>(e);
        }
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssObservationEquation *GnssReceiver::ObservationEquationList::operator()(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < eqn.size()) && (idTrans < eqn.at(idEpoch).size()))
    return eqn[idEpoch][idTrans].get();
  return nullptr;
}

/***********************************************/
/***********************************************/

void GnssReceiver::preprocessingInfo(const std::string &info, UInt countEpochs, UInt countObservations, UInt countTracks)
{
  try
  {
    if(countEpochs == NULLINDEX)
    {
      countEpochs = 0;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(useable(idEpoch))
          countEpochs++;
    }

    if(countObservations == NULLINDEX)
    {
      countObservations = 0;
      for(UInt idEpoch=0; idEpoch<observations_.size(); idEpoch++)
        for(UInt idTrans=0; idTrans<observations_.at(idEpoch).size(); idTrans++)
          if(observations_[idEpoch][idTrans])
            countObservations++;
    }

    if(countTracks == NULLINDEX)
      countTracks = tracks.size();

    preprocessingInfos.push_back(countEpochs%"%6i epochs,"s+countObservations%"%7i observations,"s+countTracks%"%4i tracks: "s+info);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::readObservations(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters,
                                    const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const Time &timeMargin, Angle elevationCutOff,
                                    const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group)
{
  try
  {
    GnssReceiverArc arc = InstrumentFile::read(fileName);
    // Counts each tracked transmitter as one observation
    preprocessingInfo("in file <"+fileName.str()+">", arc.size(),
                      std::accumulate(arc.begin(), arc.end(), UInt(0), [](UInt count, const auto &e){return count+e.satellite.size();}), 0);

    // Epochs of matched observations
    std::vector<Time> observationTimes;
    // Phase windup corrections for each transmitter
    Vector phaseWindup(transmitters.size());
    std::map<GnssType, UInt> removedTypes;

    UInt idEpoch = 0;
    // Number of deleted observations due to different reasons
    std::vector<UInt> countDeleted(4, 0);
    // loop over observation epochs
    for(UInt arcEpoch=0; arcEpoch<arc.size(); arcEpoch++)
    {
      // Search for the observation epoch corresponding to the idEpoch-th processing epoch:
      // 1) if the processing epoch falls behind the observation epoch,
      //    then this processing epoch is disabled and proceed to the next processing epoch.
      while((idEpoch < times.size()) && (times.at(idEpoch)+timeMargin < arc.at(arcEpoch).time))
        disable(idEpoch++, "missing epoch in the observation file");
      if(idEpoch >= times.size())
        break;

      // 2) if the processing epoch surpasses the observation epoch or
      //    the receiver is not usable at this epoch,
      //    then skip this observation epoch and proceed to the next observation epoch.
      if((arc.at(arcEpoch).time+timeMargin < times.at(idEpoch)) || !useable(idEpoch))
        continue;

      // 3) if the processing epoch and the observation epoch match,
      //    update the processing epoch with the observation epoch.
      times.at(idEpoch) = arc.at(arcEpoch).time;
//    clk.at(idEpoch)   = arc.at(arcEpoch).clockError;
      observationTimes.push_back(arc.at(arcEpoch).time);

      // Possibly tracked GNSS types of this receiver at this epoch
      const std::vector<GnssType> receiverTypes = definedTypes(times.at(idEpoch));

      // Counter for number of valid observations at this epoch
      UInt idObs  = 0;
      // loop over each satellite
      for(UInt k=0; k<arc.at(arcEpoch).satellite.size(); k++)
      {
        GnssType satType = arc.at(arcEpoch).satellite.at(k);
        // Index of the first obs type for this satellite
        UInt idType = 0;
        while(arc.at(arcEpoch).obsType.at(idType) != satType)
          idType++;

        // Index of this satellite in the transmitters list
        const UInt idTrans = std::distance(transmitters.begin(), std::find_if(transmitters.begin(), transmitters.end(),
                                                                              [&](auto t) {return t->PRN() == satType;}));
        // Possibly emitted GNSS types of this transmitter at this epoch
        std::vector<GnssType> transmitterTypes;
        if(idTrans < transmitters.size())
          transmitterTypes = transmitters.at(idTrans)->definedTypes(times.at(idEpoch));

        // repair GLONASS frequency number
        if((satType == GnssType::GLONASS) && transmitterTypes.size() && (transmitterTypes.front().frequencyNumber() != 9999))
          satType.setFrequencyNumber(transmitterTypes.front().frequencyNumber());

        // Valid observations from this satellite at this epoch
        GnssObservation *obs = new GnssObservation();
        // loop over each observation type for this satellite
        for(; (idType<arc.at(arcEpoch).obsType.size()) && (arc.at(arcEpoch).obsType.at(idType)==satType); idType++, idObs++)
          if((idTrans < transmitters.size()) && arc.at(arcEpoch).observation.at(idObs) && !std::isnan(arc.at(arcEpoch).observation.at(idObs)))
          {
            GnssType type = arc.at(arcEpoch).obsType.at(idType) + satType;
            // Remove frequency number for GLONASS non-G1 and non-G2 observations
            if((type == GnssType::GLONASS) && !((type == GnssType::G1) || (type == GnssType::G2)))
              type.setFrequencyNumber(9999);

            // Check the completeness of the GNSS observation type
            if(type.hasWildcard(GnssType::TYPE + GnssType::FREQUENCY + GnssType::SYSTEM + GnssType::PRN))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<" is not complete"<<Log::endl;
              continue;
            }
            if((type.frequencyNumber() == 9999) && (type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2)))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<": GLONASS frequency number not set"<<Log::endl;
              continue;
            }
            if((type.frequencyNumber() != 9999) && !((type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2))))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<": GLONASS frequency number is set"<<Log::endl;
              continue;
            }

            // Check against useType and ignoreType
            Bool use = (useType.size()==0) ? TRUE : FALSE;
            if(type.isInList(useType))
              use = TRUE;
            if(type.isInList(ignoreType))
              use = FALSE;

            // Check against receiver signal types
            if(use && receiverTypes.size() && !type.isInList(receiverTypes))
            {
              use = FALSE;
              removedTypes[type]++;
            }
            // Check against transmitter signal types
            if(use && transmitterTypes.size())
              for(const GnssType &typeTrans : GnssType::replaceCompositeSignals({type}))
                if(!typeTrans.isInList(transmitterTypes))
                {
                  use = FALSE;
                  removedTypes[type]++;
                  // Should we break here?
                }

            if(use)
              obs->push_back(GnssSingleObservation(type, arc.at(arcEpoch).observation.at(idObs)));
          }

        if(obs->size() == 0)
        {
          countDeleted[0]++;
          delete obs;
          continue;
        }
        if(idTrans >= transmitters.size())
        {
          countDeleted[1]++;
          delete obs;
          continue;
        }
        if(!obs->init(*this, *transmitters.at(idTrans), rotationCrf2Trf, idEpoch, elevationCutOff, phaseWindup(idTrans)))
        {
          countDeleted[2]++;
          delete obs;
          continue;
        }
        std::vector<GnssType> types;
        if(!obs->observationList(group, types))
        {
          countDeleted[3]++;
          delete obs;
          continue;
        }

        // store the observation in the receiver's observation list
        if(observations_.size() <= idEpoch)
          observations_.resize(idEpoch+1);
        if(observations_.at(idEpoch).size() <= idTrans)
          observations_.at(idEpoch).resize(idTrans+1, nullptr);
        if(observations_[idEpoch][idTrans])
          logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": observation already exists"<<Log::endl;
        obs->sort();
        std::swap(observations_[idEpoch][idTrans], obs);
        delete obs;
      }

      if((observations_.size() <= idEpoch) || (observations_[idEpoch].size() == 0))
        disable(idEpoch, "no useable observations found (elevationCutOff, use/ignoreTypes, defined receiver/transmitter types, missing antenna patterns)");
      idEpoch++;
    }

    for(UInt idEpoch=observations_.size(); idEpoch<times.size(); idEpoch++)
      disable(idEpoch, "missing epochs in file");

    if(removedTypes.size())
    {
      std::stringstream ss;
      for(const auto &type : removedTypes)
        ss<<"  "<<type.first.str()<<"="<<type.second;
      logWarning<<name()<<": removed undefined observations"<<ss.str()<<Log::endl;
    }

    observationSampling = medianSampling(observationTimes).seconds();
    copyObservations2ContinuousMemoryBlock();

    std::stringstream ss;
    ss << "readObservations(), deleted";
    for(UInt i=0; i<countDeleted.size(); i++)
      ss << " " << countDeleted[i];
    preprocessingInfo(ss.str());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::simulateZeroObservations(const std::vector<GnssType> &types,
                                            const std::vector<GnssTransmitterPtr> &transmitters,
                                            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, Angle elevationCutOff,
                                            const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group)
{
  try
  {
    // Phase windup corrections for each transmitter
    Vector phaseWindup(transmitters.size());
    // Number of deleted observations due to different reasons
    std::vector<UInt> countDeleted(3, 0);
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      const std::vector<GnssType> receiverTypes = definedTypes(times.at(idEpoch));

      // create observation class for each satellite
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      {
        // PRN of this satellite
        GnssType satType = transmitters.at(idTrans)->PRN();
        // list of broadcasted observation types of this satellite at this epoch
        std::vector<GnssType> transmitterTypes = transmitters.at(idTrans)->definedTypes(times.at(idEpoch));

        // set the frequency number for GLONASS satellites
        if((satType == GnssType::GLONASS) && transmitterTypes.size() && (transmitterTypes.front().frequencyNumber() != 9999))
          satType.setFrequencyNumber(transmitterTypes.front().frequencyNumber());

        GnssObservation *obs = new GnssObservation();
        // loop over each to-be-simulated observation type for this satellite
        for(UInt idType=0; idType<types.size(); idType++)
          if(types.at(idType) == satType)
          {
            GnssType type = types.at(idType) + satType;
            // remove GLONASS frequency number
            if((type == GnssType::GLONASS) && !((type == GnssType::G1) || (type == GnssType::G2)))
              type.setFrequencyNumber(9999);

            // check completeness
            if(type.hasWildcard(GnssType::TYPE + GnssType::FREQUENCY + GnssType::SYSTEM + GnssType::PRN))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<" is not complete"<<Log::endl;
              continue;
            }
            if((type.frequencyNumber() == 9999) && (type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2)))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<": GLONASS frequency number not set"<<Log::endl;
              continue;
            }
            if((type.frequencyNumber() != 9999) && !((type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2))))
            {
              logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": "<<type.str()<<": GLONASS frequency number is set"<<Log::endl;
              continue;
            }

            // check useType and ignoreType
            Bool use = (useType.size()==0) ? TRUE : FALSE;
            if(type.isInList(useType))
              use = TRUE;
            if(type.isInList(ignoreType))
              use = FALSE;

            // check against receiver signal types
            if(use && receiverTypes.size() && !type.isInList(receiverTypes))
              use = FALSE;
            // check against basic transmitter signal types
            if(use && transmitterTypes.size())
            {
              GnssType typeComposed = type;
              if(typeComposed == GnssType::PHASE) // phase can only be measured, if RANGE with is available
                typeComposed = (typeComposed & ~GnssType::TYPE) + GnssType::RANGE; // replace PHASE by RANGE (to preserve ATTRIBUTE)
              for(const GnssType &typeTrans : GnssType::replaceCompositeSignals({typeComposed}))
                if(!typeTrans.isInList(transmitterTypes))
                  use = FALSE;
            }
            if(!use)
              continue;
            // add the observation with zero value
            obs->push_back(GnssSingleObservation(type, 0.0));
          }

        if(obs->size() == 0)
        {
          countDeleted[0]++;
          delete obs;
          continue;
        }
        if(!obs->init(*this, *transmitters.at(idTrans), rotationCrf2Trf, idEpoch, elevationCutOff, phaseWindup(idTrans)))
        {
          countDeleted[1]++;
          delete obs;
          continue;
        }
        std::vector<GnssType> types;
        if(!obs->observationList(group, types))
        {
          countDeleted[2]++;
          delete obs;
          continue;
        }

        if(observations_.size() <= idEpoch)
          observations_.resize(idEpoch+1);
        if(observations_.at(idEpoch).size() <= idTrans)
          observations_.at(idEpoch).resize(idTrans+1, nullptr);
        if(observations_[idEpoch][idTrans])
          logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": observation already exists"<<Log::endl;
        obs->sort();
        std::swap(observations_[idEpoch][idTrans], obs);
        delete obs;
      }

      if((observations_.size() <= idEpoch) || (observations_[idEpoch].size() == 0))
        disable(idEpoch, "no observations simulated (elevationCutOff, use/ignoreTypes, defined receiver/transmitter types, missing antenna patterns)");
    }

    observationSampling = medianSampling(times).seconds();
    copyObservations2ContinuousMemoryBlock();

    std::stringstream ss;
    ss << "simulateObservations(), deleted";
    for(UInt i=0; i<countDeleted.size(); i++)
      ss << " " << countDeleted[i];
    preprocessingInfo(ss.str());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::simulateObservations(NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                                        const std::vector<GnssTransmitterPtr> &transmitters,
                                        const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                        const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                        UInt minObsCountPerTrack, Angle elevationTrackMinimum, GnssObservation::Group group)
{
  try
  {
    // Simulate clock error
    // --------------------
    Vector clock = noiseClock->noise(times.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      if(useable(idEpoch))
        updateClockError(idEpoch, clock(idEpoch)/LIGHT_VELOCITY);

    /** @brief A local class for a GNSS ambiguity. */
    class Ambiguity : public GnssAmbiguity
    {
    public:
      /// The GNSS types of the ambiguities
      std::vector<GnssType> types;
      /// ambiguity values, in meter
      Vector                value;

      explicit Ambiguity(GnssTrack *track, const Vector &value) : GnssAmbiguity(track), types(track->types), value(value) {}

      /** @brief Returns the ambiguity values of the given GNSS types. */
      Vector ambiguities(const std::vector<GnssType> &types) const override
      {
        Vector value(types.size());
        UInt idx;
        for(UInt idType=0; idType<types.size(); idType++)
          if(types.at(idType).isInList(this->types, idx))
            value(idType) = this->value(idx);
        return value;
      }
    };

    // init random phase ambiguities
    std::random_device randomDevice;
    // pseudo-random generator for ambiguities
    std::mt19937_64 generator;
    generator.seed(randomDevice());
    // uniform distribution for ambiguities between -10000 and 10000 cycles
    auto ambiguityRandom = std::uniform_int_distribution<Int>(-10000, 10000);

    createTracks(transmitters, minObsCountPerTrack, {});
    // simulate ambiguities for each track
    for(auto &track : tracks)
    {
      Vector value(track->types.size());
      for(UInt i=0; i<value.size(); i++)
        if(track->types.at(i) == GnssType::PHASE)
          value(i) = wavelengthFactor*track->types.at(i).wavelength() * ambiguityRandom(generator); // cycles to meter
      new Ambiguity(track.get(), value);
    }

    // reduced observations
    // --------------------
    ObservationEquationList eqnList(*this, transmitters, rotationCrf2Trf, reduceModels, group);
    removeLowElevationTracks(eqnList, elevationTrackMinimum);

    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
    {
      std::vector<GnssType> typesTrans;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        auto obs = observation(idTrans, idEpoch);
        if(obs)
          for(UInt idType=0; idType<obs->size(); idType++)
            if(!obs->at(idType).type.isInList(typesTrans))
              typesTrans.push_back(obs->at(idType).type);
      }
      if(!typesTrans.size())
        continue;

      // obs noises for each epoch and each observation type of this transmitter
      const Matrix eps = noiseObs->noise(times.size(), typesTrans.size());
      UInt idx;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(observation(idTrans, idEpoch))
        {
          const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
          GnssObservation *obs = observation(idTrans, idEpoch);
          for(UInt idType=0; idType<obs->size(); idType++)
            if(obs->at(idType).type.isInList(eqn.types, idx))
              obs->at(idType).observation = -eqn.l(idx) + eqn.sigma(idx) * eps(idEpoch, GnssType::index(typesTrans, obs->at(idType).type));
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Vector3d> GnssReceiver::estimateInitialClockErrorFromCodeObservations(const std::vector<GnssTransmitterPtr> &transmitters,
                                                                                  const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                                                                  const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                                                  Double huber, Double huberPower, Bool estimateKinematicPosition)
{
  try
  {
    // count available GNSS systems and RANGE observation types
    std::vector<GnssType> systems, types;
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        if(useable(idEpoch) && observation(idTrans, idEpoch) && observation(idTrans, idEpoch)->observationList(GnssObservation::RANGE, types))
        {
          if(!types.at(0).isInList(systems))
            systems.push_back(types.at(0) & GnssType::SYSTEM);
          break;
        }

    std::vector<Vector3d> posOld = pos;

    // number of static parameters: 3 for position (if static) + inter-system clock biases (per system)
    const UInt countStaticParameters = (estimateKinematicPosition ? 0 : 3) + systems.size()-1;
    // number of epoch parameters: 3 for position (if kinematic) + clock error
    const UInt countEpochParameters  = (estimateKinematicPosition ? 3 : 0) + 1;
    // Count the number of deleted epochs
    UInt countDisabled = 0;
    for(UInt iter=0; iter<10; iter++)
    {
      // setup observation equations: position, clock
      // --------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<Matrix> listlFull, listAFull, listBFull;
      std::vector<UInt>   listEpoch;
      std::vector<std::vector<UInt>> listIndexFull;
      // maximum number of satellites tracked at any epoch
      UInt maxSat = 0;
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        if(useable(idEpoch))
        {
          // number of observations at this epoch
          UInt                                 obsCount = 0;
          std::vector<GnssObservationEquation> eqnList;
          // count observations and setup observation equations for each transmitter
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          {
            const GnssObservation *obs = observation(idTrans, idEpoch);
            std::vector<GnssType> types;
            if(obs && obs->observationList(GnssObservation::RANGE, types))
            {
              eqnList.emplace_back(*obs, *this, *transmitters.at(idTrans), rotationCrf2Trf, reduceModels, idEpoch, TRUE, types);
              eqnList.back().eliminateGroupParameters(); // eliminate STEC
              obsCount += eqnList.back().l.rows();
            }
          }

          if(!obsCount || (eqnList.size() <= countEpochParameters))
          {
            countDisabled++;
            disable(idEpoch, "not enough observations to estimate clock errors");
            continue;
          }

          // setup combined observation equations
          Vector l(obsCount);
          Matrix A(obsCount, countStaticParameters);
          Matrix B(obsCount, countEpochParameters);

          listIndexFull.push_back({0});
          UInt idx = 0;
          for(const auto &eqn : eqnList)
          {
            const UInt count = eqn.l.rows();
            copy(eqn.l, l.row(idx, count));
            copy(eqn.A.column(GnssObservationEquation::idxClockRecv),  B.slice(idx, 0, count, 1)); // clock
            MatrixSlice Apos(estimateKinematicPosition ? B.slice(idx, 1, count, 3) : A.slice(idx, 0, count, 3));
            if(isEarthFixed())
              matMult(1., eqn.A.column(GnssObservationEquation::idxPosRecv, 3), rotationCrf2Trf(eqn.timeRecv).matrix().trans(), Apos);
            else
              copy(eqn.A.column(GnssObservationEquation::idxPosRecv, 3), Apos);
            // intersystem clock shift
            const UInt idxSys = GnssType::index(systems, eqn.types.front());
            if(idxSys > 0)
              copy(eqn.A.column(GnssObservationEquation::idxClockTrans),  A.slice(idx, estimateKinematicPosition ? (idxSys-1) : (3+idxSys-1), count, 1)); // clock bias
            idx += count;
            listIndexFull.back().push_back(idx);
          }

          listEpoch.push_back(idEpoch);
          listlFull.push_back(l);
          listBFull.push_back(B);
          maxSat = std::max(maxSat, eqnList.size());

          if(A.size())
          {
            listAFull.push_back(A);
            eliminationParameter(B, {A, l});
            listl.push_back(l);
            listA.push_back(A);
          }
        }

      if(!listEpoch.size() || (maxSat < countStaticParameters+countEpochParameters))
      {
       disable("only maximum "+maxSat%"%i satellites tracked for any epoch"s);
       return posOld;
      }

      // estimate static parameters
      // --------------------------
      // Position change if estimated statically or maximum position change of all epochs
      // if estimated kinematically in distance after each iteration
      Double maxPosDiff = 0;
      if(countStaticParameters)
      {
        // copy equations in one system
        const UInt count = std::accumulate(listl.begin(), listl.end(), UInt(0), [](UInt sum, const auto &x) {return sum+x.size();});
        Vector l(count);
        Matrix A(count, countStaticParameters);
        std::vector<UInt> index({0});
        Vector sigma;
        for(UInt i=0; i<listl.size(); i++)
        {
          const UInt idx = index.back();
          copy(listl.at(i), l.row(idx, listl.at(i).rows()));
          copy(listA.at(i), A.row(idx, listA.at(i).rows()));
          index.push_back(idx+listA.at(i).rows());
        }
        const Vector dx = Vce::robustLeastSquares(A, l, index, huber, huberPower, 30, sigma);
        for(UInt i=0; i<listEpoch.size(); i++)               // update with static parameters
          matMult(-1, listAFull.at(i), dx, listlFull.at(i));
        if(!estimateKinematicPosition)                       // udpate static position
        {
          const Vector3d dpos(dx.row(0, 3));
          maxPosDiff = dpos.r();
          for(UInt i=0; i<listEpoch.size(); i++)
            pos.at(listEpoch.at(i)) += dpos;
        }
      }

      // reconstruct epoch parameters
      // ----------------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        Vector sigma;
        const Vector y = Vce::robustLeastSquares(listBFull.at(i), listlFull.at(i), listIndexFull.at(i), huber, huberPower, 10, sigma);
        matMult(-1, listBFull.at(i), y, listlFull.at(i)); // compute residuals
        updateClockError(listEpoch.at(i), y(0)/LIGHT_VELOCITY);
        if(estimateKinematicPosition)
        {
          pos.at(listEpoch.at(i)) += Vector3d(y.row(1, 3));
          maxPosDiff = std::max(maxPosDiff, norm(y.row(1, 3)));
        }
      }

      // check convergence
      // -----------------
      if(maxPosDiff < 5.0)
        break;
    }

    preprocessingInfo("estimateInitialClockErrorFromCodeObservations(), disabled "+countDisabled%"%i epochs due to insufficient observations"s);

    // restore apriori positions and return new positions
    std::swap(pos, posOld);
    return posOld;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqnList, Double threshold, Double outlierRatio)
{
  try
  {
    // Count the number of disabled epochs
    UInt countDisabled = 0;
    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      if(useable(idEpoch))
      {
        // number of outlier satellites at this epoch
        UInt outlierCount = 0;
        // number of observed satellites at this epoch
        UInt count   = 0;
        for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          if(observation(idTrans, idEpoch))
          {
            const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
            for(UInt idType=0; idType<eqn.types.size(); idType++)
              if((eqn.types.at(idType) == GnssType::RANGE) && (std::fabs(eqn.l.at(idType)) >= threshold))
              {
                // delete all observations to a satellite at an epoch if they contain a gross code outlier
                deleteObservation(idTrans, idEpoch);
                outlierCount++;
                break;
              }
            count++;
          }

        // disable epoch if outlierRatio or more of the observed satellites have gross code outliers
        if(outlierCount >= outlierRatio * count)
        {
          disable(idEpoch, "too many gross code outliers, "+outlierCount%"%i out of "s+count%"%i transmitters"s);
          countDisabled++;
        }
      }

    preprocessingInfo("disableEpochsWithGrossCodeObservationOutliers(), disabled "+countDisabled%"%i epochs"s);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::createTracks(const std::vector<GnssTransmitterPtr> &transmitters, UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes)
{
  try
  {
    tracks.clear();
    // Count the number of deleted observations for each transmitter
    std::vector<UInt> countDeleted(transmitters.size(), 0);
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
    {
      UInt idEpochStart = 0;
      for(;;)
      {
        // find continuous track
        // ---------------------
        // find the first epoch with observations of this track
        while((idEpochStart < idEpochSize()) && !observation(idTrans, idEpochStart))
          idEpochStart++;

        GnssObservation *obs = observation(idTrans, idEpochStart);
        if(!obs) // at end?
          break;

        // phase and range types of the start epoch
        std::vector<GnssType> types;
        for(UInt idType= 0; idType <obs->size(); idType++)
          if((obs->at(idType).type == GnssType::PHASE) || (obs->at(idType).type == GnssType::RANGE))
            types.push_back(obs->at(idType).type);
        std::sort(types.begin(), types.end());

        // number of epochs in this track
        UInt countEpoch = 1;
        // find the last epoch of this track
        UInt idEpochEnd = idEpochStart;
        for(UInt idEpoch=idEpochStart+1; idEpoch<idEpochSize(); idEpoch++)
        {
          GnssObservation *obs = observation(idTrans, idEpoch);
          if(obs)
          {
            // Gap longer than 1.5 times the observation sampling ends the track
            if((times.at(idEpoch)-times.at(idEpochEnd)).seconds() > 1.5*observationSampling)
              break;

            // phase and range types of this epoch
            std::vector<GnssType> typesNew;
            for(UInt idType=0; idType<obs->size(); idType++)
              if((obs->at(idType).type == GnssType::PHASE) || (obs->at(idType).type == GnssType::RANGE))
                typesNew.push_back(obs->at(idType).type);
            // Different set of phase and range types ends the track
            if(!GnssType::allEqual(types, typesNew))
              break;

            idEpochEnd = idEpoch;
            countEpoch++;
          }
        }

        // list of frequencies of phase types in this track (additional to extraTypes)
        std::vector<GnssType> typeFrequencies;
        for(GnssType type : types)
          if((type == GnssType::PHASE) && !type.isInList(extraTypes) && !type.isInList(typeFrequencies))
            typeFrequencies.push_back(type & GnssType::FREQUENCY);

        // a valid track: minimum number of epochs and at least 2 phase frequencies
        if((countEpoch >= minObsCountPerTrack) && (typeFrequencies.size() >= 2))
        {
          tracks.push_back(std::make_shared<GnssTrack>(this, transmitters.at(idTrans).get(), idEpochStart, idEpochEnd, types));
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
            if(observation(idTrans, idEpoch))
              observation(idTrans, idEpoch)->track = tracks.back().get();
        }
        else
        {
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
          {
            deleteObservation(idTrans, idEpoch);
            countDeleted[idTrans]++;
          }
        }
        idEpochStart = idEpochEnd + 1;
      }
    }
    // total number of deleted observations for all transmitters
    UInt countDeletedTotal = std::accumulate(countDeleted.begin(), countDeleted.end(), UInt(0));
    preprocessingInfo("createTracks(), deleted "+countDeletedTotal%"%i observations"s);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::deleteTrack(UInt idTrack)
{
  try
  {
    for(UInt idEpoch=tracks.at(idTrack)->idEpochStart; idEpoch<=tracks.at(idTrack)->idEpochEnd; idEpoch++)
    {
      deleteObservation(tracks.at(idTrack)->transmitter->idTrans(), idEpoch);
      if(!useable())
        return;
    }
    tracks.erase(tracks.begin()+idTrack);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::deleteEmptyTracks()
{
  try
  {
    auto isEmpty = [](GnssTrackPtr t)
    {
      for(UInt idEpoch=t->idEpochStart; idEpoch<=t->idEpochEnd; idEpoch++)
        if(t->transmitter->useable(idEpoch) && t->receiver->observation(t->transmitter->idTrans(), idEpoch))
          return FALSE;
      return TRUE;
    };

    tracks.erase(std::remove_if(tracks.begin(), tracks.end(), isEmpty), tracks.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::removeLowElevationTracks(ObservationEquationList &eqnList, Angle minElevation)
{
  try
  {
    // count the number of deleted tracks
    UInt countDeleted = 0;
    for(UInt idTrack=tracks.size(); idTrack-->0;)
    {
      const UInt idTrans = tracks.at(idTrack)->transmitter->idTrans();
      Bool removeTrack = TRUE;
      for(UInt idEpoch=tracks.at(idTrack)->idEpochStart; idEpoch<=tracks.at(idTrack)->idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch) && (eqnList(idTrans, idEpoch)->elevationRecvAnt >= minElevation))
        {
          removeTrack = FALSE;
          break;
        }

      if(removeTrack)
      {
        countDeleted++;
        deleteTrack(idTrack);
      }
    }

    preprocessingInfo("removeLowElevationTracks(), deleted "+countDeleted%"%i tracks"s);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssTrackPtr GnssReceiver::splitTrack(ObservationEquationList &eqnList, GnssTrackPtr track, UInt idEpochSplit)
{
  try
  {
    const UInt idTrans = track->transmitter->idTrans();
    // the newly created track by splitting the old track
    GnssTrackPtr trackNew = std::make_shared<GnssTrack>(track->receiver, track->transmitter, idEpochSplit, track->idEpochEnd, track->types);
    tracks.push_back(trackNew);

    // shorten old track
    track->idEpochEnd = idEpochSplit-1;

    // connect observations to the new track
    for(UInt idEpoch=trackNew->idEpochStart; idEpoch<=trackNew->idEpochEnd; idEpoch++)
      if(observation(idTrans, idEpoch))
        observation(idTrans, idEpoch)->track = trackNew.get();

    // connect observation equations to the new track
    for(UInt idEpoch=trackNew->idEpochStart; idEpoch<=trackNew->idEpochEnd; idEpoch++)
      if(observation(idTrans, idEpoch))
        eqnList(idTrans, idEpoch)->track = trackNew.get();

    return trackNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::linearCombinations(ObservationEquationList &eqnList, GnssTrackPtr track, const std::vector<GnssType> &extraTypes,
                                      std::vector<GnssType> &typesPhase, std::vector<UInt> &idEpochs, Matrix &combinations, Double &cycles2tecu) const
{
  try
  {
    // available phase observation types in this track
    typesPhase.clear();
    idEpochs.clear();
    for(UInt idEpoch=track->idEpochStart; idEpoch<=track->idEpochEnd; idEpoch++)
      if(eqnList(track->transmitter->idTrans(), idEpoch) && observation(track->transmitter->idTrans(), idEpoch))
      {
        idEpochs.push_back(idEpoch);
        for(GnssType type : eqnList(track->transmitter->idTrans(), idEpoch)->types)
          if((type == GnssType::PHASE) && !type.isInList(typesPhase) && !type.isInList(extraTypes))
            typesPhase.push_back(type);
      }

    const Matrix Bias = GnssLambda::phaseDecorrelation(typesPhase, wavelengthFactor);
    combinations = Matrix(idEpochs.size(), Bias.columns()-1);
    UInt row = 0;
    for(UInt idEpoch : idEpochs)
    {
      const GnssObservationEquation &eqn = *eqnList(track->transmitter->idTrans(), idEpoch);
      Vector l = eqn.l;
      Matrix A(l.rows(), Bias.columns()+2);
      UInt idx;
      for(UInt idType=0; idType<eqn.types.size(); idType++) // ambiguities
        if(eqn.types.at(idType).isInList(typesPhase, idx) || (eqn.types.at(idType) == GnssType::RANGE))
        {
          l(idType) = eqn.l(idType)/eqn.sigma0(idType);
          A(idType, 0) = 1.; // range
          A(idType, 1) = eqn.types.at(idType).ionosphericFactor(); // TEC
          if(idx != NULLINDEX)
            copy(Bias.row(idx), A.slice(idType, 2, 1, Bias.columns()));
          A.row(idType) *= 1./eqn.sigma0(idType);
        }

      // skip the first, inaccurate one
      copy(leastSquares(A, l).row(2+1, Bias.columns()-1).trans(), combinations.row(row++));
    }

    // determine cycle slip size in terms of TEC
    Vector l = Bias.column(0); // cycle slips can only occur in this linear combination anymore
    Matrix A(typesPhase.size(), 2, 1.); // first column range
    for(UInt idType=0; idType<typesPhase.size(); idType++)
      A(idType, 1) = typesPhase.at(idType).ionosphericFactor(); // TEC
    cycles2tecu = std::fabs(leastSquares(A, l)(1,0)); // one cycle slip in terms of TEC
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::rangeAndTec(ObservationEquationList &eqnList, UInt idTrans, const std::vector<UInt> &idEpochs,
                               const std::vector<GnssType> &typesPhase, Vector &range, Vector &tec) const
{
  try
  {
    Matrix A(typesPhase.size(), 2, 1.);
    for(UInt k=0; k<A.rows(); k++)
      A(k, 1) = typesPhase.at(k).ionosphericFactor();

    Matrix L(typesPhase.size(), idEpochs.size());
    for(UInt i=0; i<idEpochs.size(); i++)
      for(UInt k=0; k<L.rows(); k++)
        L(k, i) = eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, typesPhase.at(k)));

    const Matrix x = leastSquares(A, L);
    range = x.row(0).trans();
    tec   = x.row(1).trans();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Double computeBias(const Vector &data, Double maxRange)
{
  try
  {
    std::vector<Double> x = data;
    std::sort(x.begin(), x.end());

    // find max length interval with data within maxRange
    auto iStartMax = x.begin();
    auto iEndMax   = x.begin();
    auto iEnd      = x.begin();
    for(auto iStart=x.begin(); iStart!=x.end(); iStart++)
    {
      while((iEnd!=x.end()) && ((*iEnd)-(*iStart) < maxRange))
        iEnd++;
      if(std::distance(iStart, iEnd) > std::distance(iStartMax, iEndMax))
      {
        iStartMax = iStart;
        iEndMax   = iEnd;
      }
    }

    //  mean of this interval
    return std::accumulate(iStartMax, iEndMax, 0.)/std::distance(iStartMax, iEndMax);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::writeTracks(const FileName &fileName, ObservationEquationList &eqnList, const std::vector<GnssType> &extraTypes) const
{
  try
  {
    if(fileName.empty())
      return;

    for(const auto &track : tracks)
      if(track->countObservations())
      {
        std::vector<GnssType> typesPhase;
        std::vector<UInt>     idEpochs;
        Matrix                combinations;
        Double                cycles2tecu;
        Vector                range, tec;
        linearCombinations(eqnList, track, extraTypes, typesPhase, idEpochs, combinations, cycles2tecu);
        rangeAndTec(eqnList, track->transmitter->idTrans(), idEpochs, typesPhase, range, tec);

        Matrix A(idEpochs.size(), 2+combinations.columns());
        // 1st column: TEC in cycles
        axpy(1./cycles2tecu, tec, A.column(1));
        // 2nd column: linear combinations
        copy(combinations, A.column(2, combinations.columns()));

        std::vector<Time> timesTrack;
        for(UInt idEpoch : idEpochs)
          timesTrack.push_back(times.at(idEpoch));
        // remove the column-wise median for 2nd and subsequent columns (columns 1..n)
        for(UInt i=1; i<A.columns(); i++)
          A.column(i) -= median(A.column(i));

        std::string typeStr;
        for(GnssType type : typesPhase)
          typeStr += type.str().substr(0, 3);
        typeStr = String::replaceAll(typeStr, "?", "");
        VariableList varList;
        varList.setVariable("station",        name());
        varList.setVariable("prn",            track->transmitter->name());
        varList.setVariable("trackTimeStart", timesTrack.front().mjd());
        varList.setVariable("trackTimeEnd",   timesTrack.back().mjd());
        varList.setVariable("types",          typeStr);
        InstrumentFile::write(fileName(varList), Arc(timesTrack, A));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsDetection(ObservationEquationList &eqnList, UInt minObsCountPerTrack, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes)
{
  try
  {
    for(UInt idTrack=0; idTrack<tracks.size(); idTrack++)
    {
      if(tracks.at(idTrack)->countObservations() >= std::max(minObsCountPerTrack, windowSize))
        cycleSlipsDetection(eqnList, tracks.at(idTrack), lambda, windowSize, tecSigmaFactor, extraTypes);
      if(tracks.at(idTrack)->countObservations() < minObsCountPerTrack)
        deleteTrack(idTrack--);
    }

    preprocessingInfo("cycleSlipsDetection()");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsDetection(ObservationEquationList &eqnList, GnssTrackPtr track, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes)
{
  try
  {
    // determine Melbourne-Wuebbena-like linear combinations
    // -----------------------------------------------------
    // list of phase types in this track (additional to extraTypes)
    std::vector<GnssType> typesPhase;
    // list of epochs in this track with valid observations
    std::vector<UInt>     idEpochs;
    // the computed linear combinations of the observations in this track
    Matrix                combinations;
    // conversion factor from cycles to TEC
    Double                cycles2tecu;
    linearCombinations(eqnList, track, extraTypes, typesPhase, idEpochs, combinations, cycles2tecu);
    // marking cycle slips at epochs
    Vector slips(idEpochs.size());
    for(UInt k=0; k<combinations.columns(); k++)
    {
      const Vector smoothed = totalVariationDenoising(combinations.column(k), lambda);
      const Double bias     = computeBias(smoothed, 0.01);
      // cycle slip if the rounded epoch-to-epoch difference of two consecutive unbiased denoised values exceeds 3/4 cycle
      for(UInt i=1; i<idEpochs.size(); i++)
        if(std::fabs(std::round(smoothed(i)-bias) - std::round(smoothed(i-1)-bias)) > 0.75)
          slips(i) = TRUE;
    }

    for(UInt i=slips.rows(); i-->0;)
      if(slips(i))
      {
        splitTrack(eqnList, track, idEpochs.at(i));
        // the original track is shortened to the epoch before the splitting epoch
        idEpochs.resize(i);
      }

    // find cycle slips in TEC based on moving window over autoregressive model residuals
    // ----------------------------------------------------------------------------------
    Vector range, tec;
    rangeAndTec(eqnList, track->transmitter->idTrans(), idEpochs, typesPhase, range, tec);

    const UInt order = 3; // AR model order
    if(windowSize && (tec.size() >= order+windowSize+1))
    {
      // Estimate AR process with Burg
      // Code mostly from books such as TimeSeriesAnalysis by James Hamilton,
      // Introduction to TimeSeries and Forecasting by brockwell and Davis,
      // and statsmodels implementaiton which seems to be the easiest way todo imo by using burg->lev

      // Compute first diff and convert to cycles
      // First diff is used to get rid of any additional issues within the tec
      // also to make sure the timeseries is stationary as possible for AR estimation
      // basically an ARIMA(3,1,0) model is used for the jump detection
      Vector x(tec.rows()-1, 0.);
      for(UInt i=1; i<tec.rows(); i++)
        x(i-1) = (tec(i) - tec(i-1))/cycles2tecu;
      x -= mean(x);

      // Burg algorithm to compute the partial autocorrelations
      Vector d(order+1, 0.);
      d(0) = 2 * quadsum(x);
      Vector pacf(order+1, 0.);
      Vector u (x.rows());
      Vector v (x.rows());

      for(UInt i=1; i<=x.rows(); i++)
      {
        u(i-1) = x(x.rows()-i);
        v(i-1) = x(x.rows()-i);
      }

      d(1) = inner(u.slice(0, u.rows()-1), u.slice(0, u.rows()-1)) + inner(v.slice(1, v.rows()-1), v.slice(1, v.rows()-1));
      pacf(1) = 2./d(1) * inner(v.slice(1, v.rows()-1), u.slice(0, u.rows()-1));

      Vector last_u(u.rows());
      Vector last_v(u.rows());
      for(UInt i=1; i<order; i++)
      {
        swap(u, last_u);
        swap(v, last_v);
        copy(last_u.slice(0, last_u.rows()-1) - pacf(i) * last_v.slice(1, last_v.rows()-1), u.slice(1, u.rows()-1));
        copy(last_v.slice(1, last_v.rows()-1) - pacf(i) * last_u.slice(0, last_u.rows()-1), v.slice(1, v.rows()-1));
        d(i+1)    = (1 - pacf(i)*pacf(i)) * d(i) - v(i)*v(i) - u(u.rows()-1)*u(u.rows()-1);
        pacf(i+1) = 2/d(i+1) * inner(v.slice(i+1, v.rows()-i-1), u.slice(i, u.rows()-i-1));
      }
      // Solve coefficients with levinson
      // first coefficient of pacf is always 1 and not necessary
      pacf = pacf.slice(1, pacf.rows()-1);

      Vector arCoeffs = pacf;
      for(UInt i=1; i<order; i++)
      {
        Vector prev = arCoeffs.slice(0, arCoeffs.rows() - (order-i));
        std::vector<UInt> prevIdx(prev.rows());
        std::iota(prevIdx.begin(), prevIdx.end(), 0);
        std::reverse(prevIdx.begin(), prevIdx.end());
        Vector prevReverse = reorder(prev, prevIdx);
        Vector temp = prev;
        temp -= arCoeffs(i) * prevReverse;
        copy(temp, arCoeffs.slice(0, arCoeffs.rows()-(order-i)));
      }

      // compute forward/backwards prediction error
      // could also be done with slicing stuff but I think this is
      // easier to understand
      // Formula: x(t) - phi_1*x(t-1) - phi_2*x(t-2) - phi_3*x(t-3) = e_forward
      //          x(t-3) - phi_1*x(t-2) - phi_2*x(t-1) - phi_3*x(t) = e_backward
      // AR coefficients should be the same according to a lot of conditions which
      // would theoretically need to be checked but honestly its not necessary
      // for this scenarios.
      Vector eForward(x.rows()-order);
      Vector eBackward(x.rows()-order);
      for(UInt i=0; i<x.rows()-order; i++)
      {
        eForward(i)  = x(i+order);
        eBackward(i) = x(x.rows()-order-i-1);
        for(UInt k=1; k<=order; k++)
        {
          eForward(i)  += (-1. * arCoeffs(k-1) * x(i+order-k));
          eBackward(i) += (-1. * arCoeffs(k-1) * x(x.rows()-i-order-1+k));
        }
      }

      std::vector<UInt> reverseIdx(eBackward.rows());
      std::iota(reverseIdx.begin(), reverseIdx.end(), 0);
      std::reverse(reverseIdx.begin(), reverseIdx.end());
      eBackward = reorder(eBackward, reverseIdx);

      // peak deteciton with the forward/backward prediction error.
      // automatic threshold scaling via MAD is used to prevent excessive
      // splitting during periods with high ionospheric variations/scintillations
      // Only if a jump in forward and backward occur a jump is decided to be true
      // in case of the first n-th windowsize samples only the backwards error decides
      // in case of the last n-th windowsize samples on the forward error decides
      // this is due to the median and MAD otherwise loosing p-order values for the MAD etc
      // Also this enables to check all values for a cycleslip except the first and last value in the ts
      // slipsdetect contain the detection value for each epoch of the obs.
      // incase a jump in forward or backward is detected will add +1 ont hat respective idx
      // Jump is concluded to be a real jump incase forward nad backward detect a jump and the reuslt is 3
      std::vector<UInt> slipsDetect(eForward.rows() + order + 1, 0);
      for(UInt i=0; i<eForward.size(); i++)
      {
        MatrixSlice currentWindowForward  = eForward.row (std::min(std::max(i, windowSize/2)-windowSize/2, eForward.rows()-windowSize), windowSize);
        MatrixSlice currentWindowBackward = eBackward.row(std::min(std::max(i, windowSize/2)-windowSize/2, eForward.rows()-windowSize), windowSize);
        const Double medForw = median(currentWindowForward);
        const Double medBack = median(currentWindowBackward);
        const Double MADForw = tecSigmaFactor*1.4826*medianAbsoluteDeviation(currentWindowForward);
        const Double MADBack = tecSigmaFactor*1.4826*medianAbsoluteDeviation(currentWindowBackward);
        const UInt   idxForw = i+order+1; // These decribe the idx of the undifferenced timeseries btw
        const UInt   idxBack = i;

        // 2 different ifs required since the forward idx is not the same as the backwards idx
        if((std::abs(eForward(i)-medForw) > 0.9) && (std::abs(eForward(i)-medForw) > MADForw))
        {
          // In case we are in the last window size only use forward prediction error as slip detection
          // doing this since the MAD is worse estimated for backwards prediction error due to missing values
          if(idxForw > slips.size()-windowSize)
            slipsDetect.at(idxForw) += 3;

          if((idxForw >= windowSize) && (idxForw <= slips.size()-windowSize))
            slipsDetect.at(idxForw) += 2;
        }
        // Since this goes backwards add +1 to the idxBack to be consistent with the jump idx from forward
        if((std::abs(eBackward(i)-medBack) > 0.9) && (std::abs(eBackward(i)-medBack) > MADBack))
        {
          // in case we are in the first window size only use backward prediction error as slip detection
          if(idxBack < windowSize)
            slipsDetect.at(idxBack+1) += 3;

          if((idxBack >= windowSize) && (idxBack <= slips.size()-windowSize))
            slipsDetect.at(idxBack+1) += 1;
        }
      }

      // list of epochs where a cycle slip is detected
      std::vector<UInt> slips;
      for(UInt i = 0; i < slipsDetect.size(); i++)
        if(slipsDetect.at(i) == 3)
          slips.push_back(i);

      for(UInt i=slips.size(); i-->0;)
      {
        splitTrack(eqnList, track, idEpochs.at(slips.at(i)));
        // shorten the original track to the epoch before the splitting epoch
        idEpochs.resize(slips.at(i));
      }
    }

    // repair cycle slips in the extra phase types (e.g. GPS L5)
    // -------------------------
    for(GnssType type : extraTypes)
    {
      const UInt   idTrans    = track->transmitter->idTrans();
      const Double wavelength = wavelengthFactor * type.wavelength();
      const Double TEC        = type.ionosphericFactor();
      for(UInt idType=0; idType<track->types.size(); idType++)
        if(track->types.at(idType) == type)
        {
          // reduce l by estimated range and tec
          Vector l(idEpochs.size());
          for(UInt i=0; i<idEpochs.size(); i++)
            l(i) = eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, track->types.at(idType))) - range(i) - TEC * tec(i);

          // fix jumps
          Double jump = 0;
          for(UInt i=1; i<l.rows(); i++)
          {
            jump += wavelength * std::round((l(i)-l(i-1))/wavelength);
            eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, track->types.at(idType))) -= jump;
            observation(idTrans, idEpochs.at(i))->at(track->types.at(idType)).observation -= jump;
          }
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList)
{
  try
  {
    // get all phase types
    std::vector<GnssType> types;
    for(const auto &track : tracks)
      for(GnssType type : track->types)
        if((type == GnssType::PHASE) && !type.isInList(types))
          types.push_back(type & ~GnssType::PRN);
    std::sort(types.begin(), types.end());

    // find two phase observations with same system and frequency
    for(UInt idType=1; idType<types.size(); idType++)
      for(UInt idType1=0; idType1<idType; idType1++)
        if(types.at(idType) == (types.at(idType1) & ~GnssType::ATTRIBUTE))
        {
          const Double wavelength = types.at(idType).wavelength();
          // compute difference
          std::vector<UInt>   idxTrans, idxEpoch;
          std::vector<Double> values;
          for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
            for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
              if(observation(idTrans, idEpoch))
              {
                const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
                UInt idx, idx1;
                if(!types.at(idType).isInList(eqn.types, idx) || !types.at(idType1).isInList(eqn.types, idx1))
                  continue;
                idxTrans.push_back(idTrans);
                idxEpoch.push_back(idEpoch);
                values.push_back((eqn.l(idx)-eqn.l(idx1))/wavelength); // diff in cycles
              }

          if(!values.size())
            continue;

          // consider bias (e.g. quarter cycles)
          Vector v0s(values.size());
          for(UInt i=0; i<values.size(); i++)
            v0s(i) = values.at(i)-std::round(values.at(i));
          const Double v0 = computeBias(v0s, 0.05);

          // fix jumps
          for(UInt i=0; i<values.size(); i++)
          {
            const Double v = wavelength * std::round(values.at(i) - v0);
            GnssObservationEquation &eqn = *eqnList(idxTrans.at(i), idxEpoch.at(i));
            eqn.l(GnssType::index(eqn.types, types.at(idType))) -= v;
            observation(idxTrans.at(i), idxEpoch.at(i))->at(types.at(idType)).observation -= v;
          }
        } // for(idType)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::trackOutlierDetection(const ObservationEquationList &eqnList, const std::vector<GnssType> &ignoreTypes, Double huber, Double huberPower)
{
  try
  {
    for(auto track : tracks)
    {
      const UInt idTrans      = track->transmitter->idTrans();
      const UInt idEpochStart = track->idEpochStart;
      const UInt idEpochEnd   = track->idEpochEnd;

      // available observations for this track
      std::vector<GnssType> types;
      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(observation(idTrans, idEpoch))
          for(const GnssType &type : eqnList(idTrans, idEpoch)->types)
            if(!type.isInList(types))
              types.push_back(type);

      // determine estimable biases (reduced by range and TEC)
      Matrix Bias = identityMatrix(types.size());
      Matrix B(types.size(), 2);
      for(UInt idType=0; idType<types.size(); idType++)
        if(types.at(idType) == GnssType::RANGE)
        {
          B(idType, 0) = 1.; //range
          B(idType, 1) = types.at(idType).ionosphericFactor(); // TEC
        }
      const Vector tau = QR_decomposition(B);
      QMult(B, tau, Bias);
      Bias = Bias.column(B.columns(), Bias.rows()-B.columns());

      // setup observation equations: range, TEC, ambiguities
      // ----------------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<UInt>   listEpoch;
      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(observation(idTrans, idEpoch))
        {
          const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);

          // observations
          Vector l = eqn.l;

          // distance and TEC
          Matrix B(l.rows(), 2);
          copy(eqn.A.column(GnssObservationEquation::idxRange), B.column(0));
          copy(eqn.A.column(GnssObservationEquation::idxSTEC),  B.column(1));

          // signal biases (includes ambiguities)
          UInt idx;
          Matrix A(l.rows(), Bias.columns());
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if(eqn.types.at(idType).isInList(types, idx))
              matMult(1., eqn.A.column(GnssObservationEquation::idxUnit+idType), Bias.row(idx), A);

          // homogenize
          for(UInt i=0; i<l.rows(); i++)
          {
            l(i)     *= 1./eqn.sigma(i);
            A.row(i) *= 1./eqn.sigma(i);
            B.row(i) *= 1./eqn.sigma(i);
          }

          // downweight ignored types
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if(eqn.types.at(idType).isInList(ignoreTypes))
            {
              l(idType)     *= 1e-3;
              A.row(idType) *= 1e-3;
              B.row(idType) *= 1e-3;
            }

          eliminationParameter(B, A, l);
          listEpoch.push_back(idEpoch);
          listl.push_back(l);
          listA.push_back(A);
        } // for(idEpoch)

      // copy equations in one system
      // ----------------------------
      const UInt count = std::accumulate(listl.begin(), listl.end(), UInt(0), [](UInt sum, const auto &x) {return sum+x.size();});
      Vector l(count);
      Matrix A(count, Bias.columns());
      std::vector<UInt> index({0});
      for(UInt i=0; i<listl.size(); i++)
      {
        const UInt idx = index.back();
        copy(listl.at(i), l.row(idx, listl.at(i).rows()));
        copy(listA.at(i), A.row(idx, listA.at(i).rows()));
        index.push_back(idx+listA.at(i).rows());
      }

      // estimate solution
      // -----------------
      Vector sigma;
      Vector x = Vce::robustLeastSquares(A, l, index, huber, huberPower, 30, sigma);

      // downweight outliers
      // -------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        const GnssObservationEquation &eqn = *eqnList(idTrans, listEpoch.at(i));
        GnssObservation &obs = *observation(idTrans, listEpoch.at(i));
        for(UInt idType=0; idType<eqn.types.size(); idType++)
          obs.at(eqn.types.at(idType)).sigma *= sigma(i);
      }

      // reduce integer ambiguities
      // --------------------------
      Vector b = Bias * x;
      for(UInt idType=0; idType<types.size(); idType++)
        if(types.at(idType) == GnssType::PHASE)
        {
          const Double lambda = types.at(idType).wavelength();
          b(idType) = lambda * std::round(b(idType)/lambda);
          for(UInt i=0; i<listEpoch.size(); i++)
          {
            eqnList(idTrans, listEpoch.at(i))->l(GnssType::index(eqnList(idTrans, listEpoch.at(i))->types, types.at(idType))) -= b(idType);
            observation(idTrans, listEpoch.at(i))->at(types.at(idType)).observation -= b(idType);
          }
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Matrix GnssReceiver::totalVariationDenoising(const_MatrixSliceRef y, Double lambda)
{
  try
  {
    Matrix x(y.rows(), y.columns());
    for(UInt col=0; col<y.columns(); col++)
    {
      // initialize total variation denoising algorithm
      UInt N = y.rows()-1;
      UInt k  = 0;                      // current sample location
      UInt k0 = 0;                      // beginning of current segment
      UInt km = 0;                      // last position where umax = -lambda
      UInt kp = 0;                      // last position where umin =  lambda
      Double vMin = y(0, col) - lambda; // lower bound for the segment's value
      Double vMax = y(0, col) + lambda; // upper bound for the segment's value
      Double uMin =  lambda;            // u is the dual variable
      Double uMax = -lambda;            // u is the dual variable

      // total variation denoising algorithm
      for(;;)
      {
        if(k == N)
        {
          x(N, col) = vMin + uMin;
          break;
        }

        if(y(k+1, col) + uMin < vMin - lambda)      // negative jump necessary
        {
          for(UInt i=k0; i<=km; i++)
            x(i, col) = vMin;
          k = k0 = km = kp = km+1;
          vMin = y(k, col);
          vMax = y(k, col) + 2*lambda;
          uMin =  lambda;
          uMax = -lambda;
        }
        else if(y(k+1, col) + uMax > vMax + lambda) // positive jump necessary
        {
          for(UInt i=k0; i<=kp; i++)
            x(i, col) = vMax;
          k = k0 = km = kp = kp+1;
          vMin = y(k, col) - 2*lambda;
          vMax = y(k, col);
          uMin =  lambda;
          uMax = -lambda;
        }
        else  // no jump necessary
        {
          k = k+1;
          uMin = uMin + y(k, col) - vMin;
          uMax = uMax + y(k, col) - vMax;
          if(uMin >= lambda)  // update of vMin
          {
            vMin = vMin + (uMin-lambda)/(k-k0+1);
            uMin = lambda;
            km = k;
          }
          if(uMax <= -lambda) // update of vMax
          {
            vMax = vMax + (uMax+lambda)/(k-k0+1);
            uMax = -lambda;
            kp = k;
          }
        }

        if(k < N)
          continue;

        if(uMin < 0.)       // vMin is too high ==> negative jump necessary
        {
          for(UInt i=k0; i<=km; i++)
            x(i, col) = vMin;
          k = k0 = km = km+1;
          vMin = y(k, col);
          uMin = lambda;
          uMax = y(k, col) + lambda - vMax;
          continue;
        }
        else if(uMax > 0.)  // vMax is too low ==> positive jump necessary
        {
          for(UInt i=k0; i<=kp; i++)
            x(i, col) = vMax;
          k = k0 = kp = kp+1;
          vMax = y(k, col);
          uMax = -lambda;
          uMin = y(k, col) - lambda - vMin;
          continue;
        }
        else
        {
          for(UInt i=k0; i<=N; i++)
            x(i, col) = vMin + uMin/(k-k0+1);
          break;
        }
      }
    }

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GnssTrack::GnssTrack(GnssReceiver *_receiver, GnssTransmitter *_transmitter, UInt _idEpochStart, UInt _idEpochEnd, const std::vector<GnssType> &_types) :
    receiver(_receiver), transmitter(_transmitter), idEpochStart(_idEpochStart), idEpochEnd(_idEpochEnd), types(_types), ambiguity(nullptr)
{
}

/***********************************************/

GnssTrack::~GnssTrack()
{
  delete ambiguity;
}

/***********************************************/

UInt GnssTrack::countObservations() const
{
  try
  {
    UInt count = 0;
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(transmitter->useable(idEpoch) && receiver->observation(transmitter->idTrans(), idEpoch))
        count++;
    return count;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssTrack::removeAmbiguitiesFromObservations(const std::vector<GnssType> &types, const std::vector<Double> &value)
{
  try
  {
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(receiver->observation(transmitter->idTrans(), idEpoch))
      {
        GnssObservation &obs = *receiver->observation(transmitter->idTrans(), idEpoch);
        UInt idx;
        for(UInt idType=0; idType<obs.size(); idType++)
          if(obs.at(idType).type.isInList(types, idx))
            obs.at(idType).observation -= value.at(idx);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
