/***********************************************/
/**
* @file parallel.h
*
* @brief Wrapper for Message Passing Interface (MPI).
* All functions are empty statements in case
* of the single processor version (parallelSingle.cpp).
*
* @author Torsten Mayer-Guerr
* @date 2004-11-13
*
*/
/***********************************************/

#ifndef __GROOPS_PARALLEL__
#define __GROOPS_PARALLEL__

#include "base/import.h"
#include "inputOutput/archiveBinary.h"
#include "inputOutput/logging.h"

/***********************************************/

class GnssType;

/***********************************************/

/** @brief Wrapper for Message Passing Interface (MPI).
* All functions are empty statements in case
* of the single processor version (parallelSingle.cpp).
* @ingroup parallelGroup */
namespace Parallel
{
  class Communicator;
  typedef std::shared_ptr<Communicator> CommunicatorPtr;

  /**
   * @brief Initializes the Communicator.
   * Must be called firstly in main.
   * @return global communicator. */
  CommunicatorPtr init(int argc, char *argv[]);

  /**
   * @brief Adds an extra communcation channel to @p comm.
   * @p receive is called on main process, if the returned send function is called by an arbitrary process.
   * This function is used for the log and must be called by every process in @a comm. */
  std::function<void(UInt type, const std::string &str)> addChannel(const std::function<void(UInt rank, UInt type, const std::string &str)> &receive, CommunicatorPtr comm);

  // =========================================================

  /**
   * @brief Splits the communicator.
   * A new group is created for each different @a color.
   * The ranks in the groups are sorted by the @a key. */
  CommunicatorPtr splitCommunicator(UInt color, UInt key, CommunicatorPtr comm);

  /**
   * @brief Creates new communicators.
   * Must be called by every process in @a comm. */
  CommunicatorPtr createCommunicator(std::vector<UInt> ranks, CommunicatorPtr comm);

  /** @brief Returns the communicator that refers to the own process only. */
  CommunicatorPtr selfCommunicator();

  // =========================================================

  /** @brief Gets the size of the communicator. */
  UInt size(CommunicatorPtr comm);

  /** @brief Gets the rank of the communicator. */
  UInt myRank(CommunicatorPtr comm);

  /** @brief Checks if the communicator is the master process (rank==0). */
  inline Bool isMaster(CommunicatorPtr comm) {return (myRank(comm) == 0);}

  /** @brief Blocks the communicator until all process have reached this routine. */
  void barrier(CommunicatorPtr comm);

  /** @brief Peeks the communicator. */
  void peek(CommunicatorPtr comm);

  /**
   * @brief Distributes exceptions thrown in @p func by a single node to all nodes.
   * Must be called by every process in @a comm.
   * Exceptions causes memory leaks due to unfinished communications.
   * Based on the idea: https://arxiv.org/abs/1804.04481 */
  void broadCastExceptions(CommunicatorPtr comm, std::function<void(CommunicatorPtr)> func);

  /** @brief Checks if the exception is caused by exceptions from other processes. */
  Bool isExternal(std::exception &e);

  // =========================================================

  /** @brief Sends raw data @a x to process with rank @a process. */
  void send(const Byte *x, UInt size, UInt process, CommunicatorPtr comm);

  /** @brief Receives raw data @a x from process with rank @a process.
  * If @a process = NULLINDEX then receive from an arbitrary process. */
  void receive(Byte *x, UInt size, UInt process, CommunicatorPtr comm);

  /** @brief Distributes raw data @a x at @a process to all other processes. */
  void broadCast(Byte *x, UInt size, UInt process, CommunicatorPtr comm);

  // =========================================================

  /** @brief Sends @a x to process with rank @a process. */
  ///@{
  template<typename T> void send(const T &x, UInt process, CommunicatorPtr comm);
  template<> void send(const UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Double   &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Time     &x, UInt process, CommunicatorPtr comm);
  template<> void send(const GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void send(const Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Receives @a x from process with rank @a process.
  * If @a process = NULLINDEX then receive from an arbitrary process. */
  ///@{
  template<typename T> void receive(T &x, UInt process, CommunicatorPtr comm);
  template<> void receive(UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Double   &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Time     &x, UInt process, CommunicatorPtr comm);
  template<> void receive(GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void receive(Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Distributes @a x at @a process to all other processes. */
  ///@{
  template<typename T> void broadCast(T &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(UInt     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Double   &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Bool     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Angle    &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Time     &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(GnssType &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Vector3d &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Vector   &x, UInt process, CommunicatorPtr comm);
  template<> void broadCast(Matrix   &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Sums up @a x at all processes (also rank 0) and sends the result to @a process. */
  ///@{
  void reduceSum(UInt                &x, UInt process, CommunicatorPtr comm);
  void reduceSum(Double              &x, UInt process, CommunicatorPtr comm);
  void reduceSum(Bool                &x, UInt process, CommunicatorPtr comm);
  void reduceSum(Matrix              &x, UInt process, CommunicatorPtr comm);
  void reduceSum(std::vector<Double> &x, UInt process, CommunicatorPtr comm);
  ///@}

  /** @brief Finds min/max of @a x at all processes (also rank 0) and sends the result to @a process. */
  ///@{
  void reduceMin(UInt   &x, UInt process, CommunicatorPtr comm);
  void reduceMin(Double &x, UInt process, CommunicatorPtr comm);
  void reduceMax(UInt   &x, UInt process, CommunicatorPtr comm);
  void reduceMax(Double &x, UInt process, CommunicatorPtr comm);
  ///@}

  // =========================================================

  /**
   * @brief Parallelizes a loop.
   * Calls the function @a func for every index in [0, @a count ).
   * The different calls are distributed among all processes except the master.
  * @return The process number for @a i is returned (valid at master). */
  template<typename T> std::vector<UInt> forEach(UInt count, T func, CommunicatorPtr comm, Bool timing=TRUE);

  /**
   * @brief Parallelizes a loop.
   * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
   * The different calls are distributed among all processes except the master.
   * The result in @a vec is only valid at master.
   * @return The process number for @a i is returned (valid at master). */
  template<typename A, typename T> std::vector<UInt> forEach(std::vector<A> &vec, T func, CommunicatorPtr comm, Bool timing=TRUE);

  /**
   * @brief Parallelizes a loop.
   * Calls @a func(i) for every @a i in [0,count).
   * The different calls are distributed among all processes except the master.
   * @return The process number for @a i is returned (valid at master). */
  template<typename T> std::vector<UInt> forEachInterval(UInt count, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing=TRUE);

  /**
   * @brief Parallelizes a loop.
   * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
   * The different calls are distributed among all processes except the master.
   * The result in @a vec is only valid at master.
   * @return The process number for @a i is returned (valid at master). */
  template<typename A, typename T> std::vector<UInt> forEachInterval(std::vector<A> &vec, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing=TRUE);

  /**
   * @brief Parallelizes a loop.
   * Calls @a func(i) for every @a i in [0,count).
   * The different calls are distributed using @a processNo (without master).
   * The result in @a vec is only valid at master. */
  template<typename T> void forEachProcess(UInt count, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing=TRUE);

  /**
   * @brief Parallelizes a loop.
   * Calls @a vec[i]=func(i) for every @a i in [0,vec.size()).
   * The different calls are distributed using @a processNo (without master).
   * The result in @a vec is only valid at master. */
  template<typename A, typename T> void forEachProcess(std::vector<A> &vec, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing=TRUE);
} // end namespace Parallel

/***********************************************/

/** @brief Loop with timing.
* @ingroup parallelGroup */
namespace Single
{
  /**
   * @brief loop with timing.
   * Calls @a func(i) for every @a i in [0,count). */
  template<typename T> void forEach(UInt count, T func, Bool timing=TRUE);
} // end namespace Single

/***********************************************/
/***** INLINES ***********************************/
/***********************************************/

template<typename T>
inline void Parallel::send(const T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  std::stringstream stream; //(std::ios::binary);
  OutArchiveBinary oa(stream, "", MAX_UINT);
  oa<<nameValue("xxx", x);
  std::string str = stream.str();
  UInt size = str.size();
  send(size, process, comm);
  send(str.data(), size, process, comm);
}

/***********************************************/

template<typename T>
inline void Parallel::receive(T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  UInt size;
  receive(size, process, comm);
  Byte *str = new Byte[size+1];
  receive(str, size, process, comm);
  std::stringstream stream(std::string(str, size)); //, std::ios::binary);
  InArchiveBinary ia(stream);
  ia>>nameValue("xxx", x);
  delete[] str;
}

/***********************************************/

template<typename T>
inline void Parallel::broadCast(T &x, UInt process, CommunicatorPtr comm)
{
  if(size(comm)<=1)
    return;
  if(Parallel::myRank(comm) == process)
  {
    std::stringstream stream(std::ios_base::out | std::ios::binary);
    OutArchiveBinary oa(stream, "", MAX_UINT);
    oa<<nameValue("xxx", x);
    std::string str = stream.str();
    UInt size = str.size();
    broadCast(size, process, comm);
    broadCast(const_cast<Byte*>(str.data()), size, process, comm);
  }
  else
  {
    UInt size;
    broadCast(size, process, comm);
    Byte *str = new Byte[size+1];
    broadCast(str, size, process, comm);
    std::stringstream stream(std::string(str, size), std::ios_base::in | std::ios::binary);
    InArchiveBinary ia(stream);
    ia>>nameValue("xxx", x);
    delete[] str;
  }
}

/***********************************************/
/***********************************************/

template<typename T>
inline std::vector<UInt> Parallel::forEach(UInt count, T func, CommunicatorPtr comm, Bool timing)
{
  return forEachInterval(count, {0, count}, func, comm, timing);
}

/***********************************************/

template<typename A, typename T>
inline std::vector<UInt> Parallel::forEach(std::vector<A> &vec, T func, CommunicatorPtr comm, Bool timing)
{
  return forEachInterval(vec, {0, vec.size()}, func, comm, timing);
}

/***********************************************/

template<typename T>
inline std::vector<UInt> Parallel::forEachInterval(UInt count, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::vector<UInt> processNo(count, 0);

    // single process version
    // ----------------------
    if(size(comm) < 3)
    {
      // single process version
      if(isMaster(comm))
      {
        Log::Timer timer(count, 1, timing);
        for(UInt i=0; i<count; i++)
        {
          timer.loopStep(i);
          func(i);
        }
        timer.loopEnd();
      }
      return processNo;
    }

    if(count!=interval.back())
      throw(Exception("interval size and count differ"));

    // parallel version
    // ----------------
    if(isMaster(comm))
    {
      std::vector<UInt> countInInterval(interval.size()-1, 0);
      std::vector<UInt> processedInterval(size(comm), NULLINDEX);

      // master distributes the loop numbers
      UInt process, index;
      Log::Timer timer(count, size(comm)-1, timing);
      for(UInt i=0; i<count; i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process

        // can we compute func in the same interval?
        UInt idInterval = processedInterval.at(process);
        if((idInterval==NULLINDEX) || (countInInterval.at(idInterval) >= interval.at(idInterval+1)-interval.at(idInterval)))
        {
          // search new interval to compute
          UInt maxLeft = 0;
          for(UInt k=0; k<countInInterval.size(); k++)
          {
            UInt left = interval.at(k+1)-interval.at(k)-countInInterval.at(k);
            // interval not used?
            if((countInInterval.at(k) == 0) && (left>0))
            {
              idInterval = k;
              break;
            }
            if(left>maxLeft)
            {
              maxLeft = left;
              idInterval = k;
            }
          }
          processedInterval.at(process) = idInterval;
        }

        UInt id = interval.at(idInterval) + countInInterval.at(idInterval);
        countInInterval.at(idInterval)++;
        send(id, process, comm);       // send new loop number to be computed at process
        processNo.at(id) = process;
        timer.loopStep(i);
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<size(comm); i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process
        send(NULLINDEX, process, comm);    // end signal
      }
      timer.loopEnd();
    }
    else // clients
    {
      send(myRank(comm), 0, comm);
      send(NULLINDEX, 0, comm); // no results computed yet
      for(;;)
      {
        UInt i;
        receive(i,0, comm);
        if(i==NULLINDEX)
          break;
        func(i);
        send(myRank(comm), 0, comm);
        send(i, 0, comm);
      }
    }

    broadCast(processNo, 0, comm);
    return processNo;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename A, typename T>
inline std::vector<UInt> Parallel::forEachInterval(std::vector<A> &vec, const std::vector<UInt> &interval, T func, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::vector<UInt> processNo(vec.size(), 0);

    // single process version
    // ----------------------
    if(size(comm) < 3)
    {
      // single process version
      if(isMaster(comm))
      {
        Log::Timer timer(vec.size(), 1, timing);
        for(UInt i=0; i<vec.size(); i++)
        {
          timer.loopStep(i);
          vec[i] = func(i);
        }
        timer.loopEnd();
      }
      return processNo;
    }

    if(vec.size()!=interval.back())
      throw(Exception("interval size and vec.size() differ"));

    // parallel version
    // ----------------
    if(isMaster(comm))
    {
      std::vector<UInt> countInInterval(interval.size()-1, 0);
      std::vector<UInt> processedInterval(size(comm), NULLINDEX);

      // master distributes the loop numbers
      UInt process, index;
      Log::Timer timer(vec.size(), size(comm)-1, timing);
      for(UInt i=0; i<vec.size(); i++)
      {
        receive(process, NULLINDEX, comm); // which process needs work?
        receive(index,   process, comm);  // loop numer be computed at process
        if(index!=NULLINDEX)
        {
          receive(vec[index], process, comm); // receive result
        }

        // can we compute func in the same interval?
        UInt idInterval = processedInterval.at(process);
        if((idInterval==NULLINDEX) || (countInInterval.at(idInterval) >= interval.at(idInterval+1)-interval.at(idInterval)))
        {
          // search new interval to compute
          UInt maxLeft = 0;
          for(UInt k=0; k<countInInterval.size(); k++)
          {
            UInt left = interval.at(k+1)-interval.at(k)-countInInterval.at(k);
            // interval not used?
            if((countInInterval.at(k) == 0) && (left>0))
            {
              idInterval = k;
              break;
            }
            if(left>maxLeft)
            {
              maxLeft = left;
              idInterval = k;
            }
          }
          processedInterval.at(process) = idInterval;
        }

        UInt id = interval.at(idInterval) + countInInterval.at(idInterval);
        countInInterval.at(idInterval)++;
        send(id, process, comm);           // send new loop number to be computed at process
        processNo.at(id) = process;
        timer.loopStep(i);
      }
      // send to all processes the end signal (NULLINDEX)
      for(UInt i=1; i<size(comm); i++)
      {
        receive(process, NULLINDEX, comm);    // which process needs work?
        receive(index,   process, comm);      // loop numer be computed at process
        if(index!=NULLINDEX)
          receive(vec[index], process, comm); // receive result
        send(NULLINDEX, process, comm);       // end signal
      }
      timer.loopEnd();
    }
    else // clients
    {
      send(myRank(comm), 0, comm);
      send(NULLINDEX, 0, comm); // no results computed yet
      for(;;)
      {
        UInt i;
        receive(i,0, comm);
        if(i==NULLINDEX)
          break;
        vec[i] = func(i);
        send(myRank(comm), 0, comm);
        send(i, 0, comm);
        send(vec[i], 0, comm);
      }
    }

    broadCast(processNo, 0, comm);
    return processNo;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename T>
inline void Parallel::forEachProcess(UInt count, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::set<UInt> procs;
    for(UInt p : processNo)
      procs.insert(p);

    UInt idx;
    Log::Timer timer(count, procs.size(), timing);
    for(UInt i=0; i<count; i++)
    {
      timer.loopStep(i);
      if(myRank(comm) == processNo.at(i))
      {
        func(i);
        if(!isMaster(comm)) send(i, 0, comm);
      } // if(arcs.at(i))
      else if(isMaster(comm))
      {
        receive(idx, NULLINDEX, comm);
      }
    } // for(i)
    barrier(comm);
    timer.loopEnd();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<typename A, typename T>
inline void Parallel::forEachProcess(std::vector<A> &vec, T func, const std::vector<UInt> &processNo, CommunicatorPtr comm, Bool timing)
{
  try
  {
    std::set<UInt> procs;
    for(UInt p : processNo)
      procs.insert(p);

    UInt idx = 0;
    Log::Timer timer(vec.size(), procs.size(), timing);
    for(UInt i=0; i<vec.size(); i++)
    {
      timer.loopStep(i);
      if(myRank(comm) == processNo.at(i))
      {
        vec[i] = func(i);
        if(!isMaster(comm)) send(i, 0, comm);
        if(!isMaster(comm)) send(vec[i], 0, comm);
      } // if(arcs.at(i))
      else if(isMaster(comm))
      {
        receive(idx, NULLINDEX, comm);
        receive(vec[idx], processNo.at(idx), comm);
      }
    } // for(i)
    barrier(comm);
    timer.loopEnd();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<typename T>
inline void Single::forEach(UInt count, T func, Bool timing)
{
{
  try
  {
    Log::Timer timer(count, 1, timing);
    for(UInt i=0; i<count; i++)
    {
      timer.loopStep(i);
      func(i);
    } // for(i)
    timer.loopEnd();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
}

/***********************************************/

#endif /* __GROOPS__ */
