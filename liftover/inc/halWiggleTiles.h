/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALWIGGLETILES_H
#define _HALWIGGLETILES_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "hal.h"

namespace hal {

/** Memory structure to keep track of wiggle results by tiling the genome
 * into regular intervals.  The idea is that if we are only writing a 
 * subregion, then we don't bother allocating space for the whole genome 
*/
template <class T>
class WiggleTiles
{
public:
   
   WiggleTiles();
   virtual ~WiggleTiles();

   /** Initialize the dimensions of the structure and default value.  No
    * memory is allocated */
   void init(hal_size_t genomeSize, T defualtValue, 
             hal_size_t tileSize);
   void clear();

   /** Get a value.  If the position does not exist in a tile, then we return
    * the default value */
   T get(hal_index_t pos) const;

   /** Set a value.  If the position does not exist in a tile, then we allocate
    * the appropriate tile.  The tile is filled with a default value except
    * pos which is set to val*/
   void set(hal_index_t pos, T val);

   /** Test if a value was written to using set().  Doing a get() isn't 
    * really sufficient because it will return the default value in the case
    * where it was not set, or if it was set with the default value */
   bool exists(hal_index_t pos) const;

   /** Methods to get basic structure info */
   hal_size_t getGenomeSize() const;
   hal_size_t getTileSize() const;
   hal_size_t getNumTiles() const;
   bool isTileEmpty(hal_size_t tile) const;
   T getDefaultValue() const;
   
protected:

   std::vector<std::vector<T> > _tiles;
   std::vector<std::vector<bool> > _bits;
   hal_size_t _tileSize;
   hal_size_t _genomeSize;
   hal_size_t _lastTileSize;
   T _defaultValue;
};


// INLINE METHODS
template <class T>
inline WiggleTiles<T>::WiggleTiles() : _tileSize(0), _genomeSize(0),
                                    _lastTileSize(0)
{
}

template <class T>
inline WiggleTiles<T>::~WiggleTiles() 
{
}

template <class T>
inline void WiggleTiles<T>::init(hal_size_t genomeSize, T defaultValue, 
                         hal_size_t tileSize)
{
  _genomeSize = genomeSize;
  _defaultValue = defaultValue;
  _tileSize = std::min(tileSize, genomeSize);
  _lastTileSize = genomeSize % _tileSize;
  hal_size_t numTiles = genomeSize / _tileSize;
  if (_lastTileSize > 0)
  {
    ++numTiles;
  }
  else
  {
    _lastTileSize = _tileSize;
  }
  for (size_t i = 0; i < _tiles.size(); ++i)
  {
    _tiles[i].clear();
    _bits[i].clear();
  }
  _tiles.resize(numTiles);
  _bits.resize(numTiles);
}

template <class T>
inline void WiggleTiles<T>::clear()
{
  _tiles.clear();
  _bits.clear();
  _tileSize = 0;
  _genomeSize = 0;
  _lastTileSize = 0;
}

template <class T>
inline T WiggleTiles<T>::get(hal_index_t pos) const
{
  assert(pos < _genomeSize);
  hal_size_t tile = pos / _tileSize;
  assert(tile < _tiles.size());
  assert(_tiles[tile].size() == 0 || _tiles[tile].size() == _tileSize ||
         _tiles[tile].size() == _lastTileSize);
  if (_tiles[tile].size() == 0)
  {
    return _defaultValue;
  }
  hal_size_t offset = pos % _tileSize;
  return _tiles[tile][offset];
}

template <class T>
inline void WiggleTiles<T>::set(hal_index_t pos, T val)
{
  assert(pos < _genomeSize);
  hal_size_t tile = pos / _tileSize;
  assert(tile < _tiles.size());
  assert(_tiles[tile].size() == 0 || _tiles[tile].size() == _tileSize ||
         _tiles[tile].size() == _lastTileSize);
  if (_tiles[tile].size() == 0)
  {
    hal_size_t len = tile == _tiles.size() - 1 ? _lastTileSize : _tileSize;
    _tiles[tile].assign(len, _defaultValue);
    _bits[tile].assign(len, false);
  }
  hal_size_t offset = pos % _tileSize;
  _tiles[tile][offset] = val;
  _bits[tile][offset] = true;
}

template <class T>
inline bool WiggleTiles<T>::exists(hal_index_t pos) const
{
  assert(pos < _genomeSize);
  hal_size_t tile = pos / _tileSize;
  assert(tile < _tiles.size());
  assert(_tiles[tile].size() == 0 || _tiles[tile].size() == _tileSize ||
         _tiles[tile].size() == _lastTileSize);
  if (_tiles[tile].size() == 0)
  {
    return false;
  }
  hal_size_t offset = pos % _tileSize;
  return _bits[tile][offset];
}

template <class T>
inline hal_size_t WiggleTiles<T>::getGenomeSize() const
{
  return _genomeSize;
}

template <class T>
inline hal_size_t WiggleTiles<T>::getTileSize() const
{
  return _tileSize;
}

template <class T>
inline hal_size_t WiggleTiles<T>::getNumTiles() const
{
  return _tiles.size();
}

template <class T>
inline bool WiggleTiles<T>::isTileEmpty(hal_size_t tile) const
{
  assert(tile < _tiles.size());
  return _tiles[tile].empty();
}

template <class T>
inline T WiggleTiles<T>::getDefaultValue() const
{
  return _defaultValue;
}

}
#endif
