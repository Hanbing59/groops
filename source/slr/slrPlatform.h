/***********************************************/
/**
* @file slrPlatform.h
*
* @brief SLR station or satellite.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPLATFORM__
#define __GROOPS_SLRPLATFORM__

#include "files/filePlatform.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrPlatform;
typedef std::shared_ptr<SlrPlatform> SlrPlatformPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for a SLR station or satellite */
class SlrPlatform
{
   /// Whether the platform is usable at a given epoch (or all epochs).
   Bool     useable_;

public:
   /// ID of the platform (set by Slr::init())
   UInt     id_;
   Platform platform;
   Double   rangeBias;

public:
  /// Constructor.
  SlrPlatform(const Platform &platform);

  /// Destructor.
  virtual ~SlrPlatform() {}

  /** @brief Gets the name of the platform. */
  std::string name() const {return platform.name;}

  /** @brief Checks if the platform is usable at a given epoch (or all epochs). */
  Bool useable() const {return useable_;}

  /** @brief Disables platform. */
  virtual void disable(const std::string &reason);

  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

inline SlrPlatform::SlrPlatform(const Platform &platform)
  : useable_(TRUE), platform(platform), rangeBias(0.) {}


/***********************************************/

inline void SlrPlatform::disable(const std::string &/*reason*/)
{
  try
  {
    useable_ = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrPlatform::save(OutArchive &oa) const
{
  oa<<nameValue("useable", useable_);
}

/***********************************************/

inline void SlrPlatform::load(InArchive  &ia)
{
  ia>>nameValue("useable", useable_);
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
