/*
Minetest Canyon Mapgen
Copyright (C) 2010-2015 kwolekr, Ryan Kwolek <kwolekr@minetest.net>
Copyright (C) 2010-2015 paramat, Matt Gregory
Copyright (C) 2016 MillersMan, Johannes Matokic

makeChunk based on other mapgens by kwolekr and paramat.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef MAPGENCANYON_H
#define MAPGENCANYON_H

#include "mapgen.h"

class BiomeGenOriginal;

struct MapgenCanyonParams : public MapgenSpecificParams
{
	u32 spflags;
	int world_level;
	// Controls how many streams there are and how big rivers can become
	// Typical range is from 0.00001 to 0.001. Higher values not recommendet.
	float rain_factor;
	// Factor of how far the rivers dominance of climate spreads.
	float river_dominance_expansion;
	// Controls how much river adjusts to its environment instead of enforcing it.
	// 0: Will crash
	// ]0,1[: Rivers will try to go through the mountain tops
	// 1: River routing is decoupled from landscape
	// > 1: Rivers try to flow through existing valleys
	float river_fitting;
	// Controls how much rivers should avoid cell-borders.
	// Rivers at cell border result in straight rivers with cliffs and almost
	// touching other rivers. (1 = off; ~2.7 = extreme)
	float river_nice_routing;
	// Controls how much rivers try to flow at diagonals and avoid sharp turns.
	// (1 = disabled; ~4 = extreme)
	// NOTE: It's not sure whether it actually improves the world
	float river_rounding;
	// Controls how much rivers will flatten the environment.
	// 1 = deep canyons. 50 = mostly flat land
	// Note: high river_valley_broadening improves river_fitting result
	float river_valley_broadening;
	// The ratio of a rivers width to depth. 
	float river_width_to_depth_ratio;
	float ground_base_scale;
	float ground_noise_scale;
	float ground_noise_off;

	float cave_width;
	NoiseParams np_filler_depth;
	NoiseParams np_cave1;
	NoiseParams np_cave2;

	MapgenCanyonParams();
	void readParams(const Settings *settings);
	void writeParams(Settings *settings) const;
};

enum {
	// Content size if the size of usable data in chunk
	CanyonChunkContentSize = 64,
	// Area at border of chunk that gets polluted at every chunk division
	CanyonChunkJunkMargin = 4,
	// Border size is for extra nodes that are needed to avoid damage of content
	// (1x for X-split, 1x for Y-split, 2x for restoring full border after split)
	CanyonChunkBorderSize = 4 * CanyonChunkJunkMargin,
	// Buffer size is the size for the whole buffer
	CanyonChunkBufferSize = CanyonChunkContentSize + 2 * CanyonChunkBorderSize
};

struct CanyonNode
{
	// Fluid-level > 0 is land.
	// Fluid-level < 0 is sea.
	float fluid_level;
	float ground_level;
	float flow_x;
	float flow_y;
	float river_dominance;
	// Parents are: 0=-x, 1=-y, 2=+x, 3=+y
	unsigned int parent : 2;
	bool join_x : 1;
	bool join_y : 1;
};

struct CanyonPatch
{
	int level;
	v2s16 pos;
	CanyonNode data[CanyonChunkBufferSize * CanyonChunkBufferSize]; 
};

class MapgenCanyon : public MapgenBasic
{
	enum { Border = CanyonChunkBorderSize };
public:
	MapgenCanyon(int mapgenid, MapgenParams* params, EmergeManager* emerge);
	~MapgenCanyon();
public:
	void makeChunk(BlockMakeData* data);
	int getSpawnLevelAtPoint(v2s16 p);
private:
	// Gets the patch that contains the given coordinate
	CanyonPatch *getPatch(int x, int y, int level);
	CanyonPatch *getChunk(v2s32 pos, int level);
	CanyonPatch *createRoot();
	CanyonPatch *divide(const CanyonPatch *parent, bool right_side);
	float groundNoise(float base, int x, int y, int level);
private:
	BiomeGenOriginal *m_bgen;

	int world_level;
	float rain_factor;
	float river_dominance_expansion;
	float river_fitting;
	float river_nice_routing;
	float river_rounding;
	float river_valley_broadening;
	float river_width_to_depth_ratio;
	float ground_base_scale;
	float ground_noise_scale;
	float ground_noise_off;
	content_t c_sand;
	
	// TODO Limit maximum number of cached chunks
	typedef std::map<v3s32, CanyonPatch*> ChunkStorage;
	ChunkStorage m_chunks;
};

class MapgenFactoryCanyon : public MapgenFactory {
public:
	Mapgen *createMapgen(int mgid, MapgenParams *params, EmergeManager *emerge)
	{
		return new MapgenCanyon(mgid, params, emerge);
	}
	MapgenSpecificParams *createMapgenParams()
	{
		return new MapgenCanyonParams();
	}
};

#endif // MAPGENCANYON_H
