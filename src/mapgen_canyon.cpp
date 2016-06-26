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

#include "mapgen_canyon.h"
#include "emerge.h"
#include "mg_biome.h"
#include "mg_decoration.h"
#include "mg_ore.h"
#include "settings.h"

#include <iostream>

static FlagDesc flagdesc_mapgen_canyon[] = {
	{NULL, 0}
};

struct CanyonBaseCursor {
	const CanyonNode *pos;
	inline CanyonBaseCursor(const CanyonNode *n)
		: pos(n)
	{ }
	inline static int startOffset(int x, int y, bool right)
	{
		const int extra = (right ? CanyonChunkContentSize / 2 : 0) +
				CanyonChunkBorderSize / 2;
		return (extra + x) * CanyonChunkBufferSize + y;
	}
	inline CanyonBaseCursor(const CanyonPatch &patch, int x, int y, bool right)
		: pos(patch.data + startOffset(x, y, right))
	{ }
	inline CanyonBaseCursor at(int dx, int dy) const
	{
		return CanyonBaseCursor((pos + dx * CanyonChunkBufferSize + dy));
	}
	inline void move(int dx, int dy) { pos += dx * CanyonChunkBufferSize + dy; }
	inline float flow_nx() const { return pos->flow_y; }
	inline float flow_ny() const { return pos->flow_x; }
	inline float flow_px() const { return -at(1, 0).flow_nx(); }
	inline float flow_py() const { return -at(0, 1).flow_ny(); }
	inline float flow_y() const { return flow_ny() - flow_py(); }
    inline float flow() const
    {
		return std::max(0.f, flow_nx()) + std::max(0.f, flow_ny()) +
				std::max(0.f, flow_px()) + std::max(0.f, flow_py());
	}
	inline bool join_nx() const { return pos->join_y; }
	inline bool join_ny() const { return pos->join_x; }
	inline bool join_px() const { return at(1, 0).join_nx(); }
	inline bool join_py() const { return at(0, 1).join_ny(); }
	inline float fluid_level() const { return pos->fluid_level; }
	inline float ground_level() const { return pos->ground_level; }
	inline float river_dominance() const { return pos->river_dominance; }
	inline int parent() const { return pos->parent ^ 0x1; }
};

struct CanyonEditCursor {
	CanyonNode *pos;
	inline CanyonEditCursor(CanyonNode *n)
		: pos(n)
	{ }
	inline CanyonEditCursor(CanyonPatch &patch, int x, int y)
		: pos(patch.data + (y * CanyonChunkBufferSize + x * 2))
	{ }
	inline CanyonEditCursor at(int dx, int dy) const
	{
		return CanyonEditCursor((pos + dy * CanyonChunkBufferSize + dx * 2));
	}
	inline void move(int dx, int dy) { pos += dy * CanyonChunkBufferSize + dx * 2; }
	inline float left_flow_nx() const { return pos->flow_x; }
	inline void set_left_flow_nx(float v) { pos->flow_x = v; }
	inline float right_flow_nx() const { return (pos + 1)->flow_x; }
	inline void set_right_flow_nx(float v) { (pos + 1)->flow_x = v; }
	inline float left_flow_ny() const { return pos->flow_y; }
	inline void set_left_flow_ny(float v) { pos->flow_y = v; }
	inline float right_flow_ny() const { return (pos + 1)->flow_y; }
	inline void set_right_flow_ny(float v) { (pos + 1)->flow_y = v; }
	inline float left_flow_px() const { return -right_flow_nx(); }
	inline void set_left_flow_px(float v) { set_right_flow_nx(-v); }
	inline float right_flow_px() const { return -at(1,0).left_flow_nx(); }
	inline void set_right_flow_px(float v) { at(1,0).set_left_flow_nx(-v); }
	inline float left_flow_py() const { return -at(0,1).left_flow_ny(); }
	inline void set_left_flow_py(float v) { at(0,1).set_right_flow_ny(-v); }
	inline float right_flow_py() const { return -at(0,1).right_flow_ny(); }
	inline void set_right_flow_py(float v) { at(0,1).set_right_flow_ny(-v); }
	inline bool left_join_nx() const { return pos->join_x; }
	inline void set_left_join_nx(bool j) const { pos->join_x = j; }
	inline bool right_join_nx() const { return (pos + 1)->join_x; }
	inline void set_right_join_nx(bool j) const { (pos + 1)->join_x = j; }
	inline bool left_join_ny() const { return pos->join_y; }
	inline void set_left_join_ny(bool j) const { pos->join_y = j; }
	inline bool right_join_ny() const { return (pos + 1)->join_y; }
	inline void set_right_join_ny(bool j) const { (pos + 1)->join_y = j; }
	inline bool left_join_px() const { return right_join_nx(); }
	inline void set_left_join_px(bool j) const { set_right_join_nx(j); }
	inline bool right_join_px() const { return at(1, 0).left_join_nx(); }
	inline void set_right_join_px(bool j) const { at(1, 0).set_left_join_nx(j); }
	inline bool left_join_py() const { return at(0, 1).left_join_ny(); }
	inline void set_left_join_py(bool j) const { at(0, 1).set_left_join_ny(j); }
	inline bool right_join_py() const { return at(0, 1).right_join_ny(); }
	inline void set_right_join_py(bool j) const { at(0, 1).set_right_join_ny(j); }
	inline float left_fluid_level() const { return pos->fluid_level; }
	inline void set_left_fluid_level(float v) { pos->fluid_level = v; }
	inline float right_fluid_level() const { return (pos + 1)->fluid_level; }
	inline void set_right_fluid_level(float v) { (pos + 1)->fluid_level = v; }
	inline float left_ground_level() const { return pos->ground_level; }
	inline void set_left_ground_level(float v) { pos->ground_level = v; }
	inline float right_ground_level() const { return (pos + 1)->ground_level; }
	inline void set_right_ground_level(float v) { (pos + 1)->ground_level = v; }
	inline float left_river_dominance() const { return pos->river_dominance; }
	inline void set_left_river_dominance(float v) { pos->river_dominance = v; }
	inline float right_river_dominance() const { return (pos + 1)->river_dominance; }
	inline void set_right_river_dominance(float v) { (pos + 1)->river_dominance = v; }
	inline int left_parent() const { return pos->parent; }
	inline void set_left_parent(int p) { pos->parent = p; }
	inline int right_parent() const { return (pos + 1)->parent; }
	inline void set_right_parent(int p) { (pos + 1)->parent = p; }
};

MapgenCanyon::MapgenCanyon(int mapgenid, MapgenParams* params, EmergeManager* emerge)
	: MapgenBasic(mapgenid, params, emerge)
{
	// NOTE: MapgenCanyon has a hard dependency on BiomeGenOriginal
	this->m_bgen = (BiomeGenOriginal *)biomegen;
	
	MapgenCanyonParams *sp = (MapgenCanyonParams*)params->sparams;
	this->spflags          = sp->spflags;
	this->world_level = sp->world_level;
	this->rain_factor = sp->rain_factor;
	this->river_dominance_expansion = sp->river_dominance_expansion;
	this->river_fitting = sp->river_fitting;
	this->river_nice_routing = sp->river_nice_routing;
	this->river_rounding = sp->river_rounding;
	this->river_valley_broadening = sp->river_valley_broadening;
	this->river_width_to_depth_ratio = sp->river_width_to_depth_ratio;
	this->ground_base_scale = sp->ground_base_scale;
	this->ground_noise_scale = sp->ground_noise_scale;
	this->ground_noise_off = sp->ground_noise_off;
	this->cave_width         = sp->cave_width;

	noise_filler_depth = new Noise(&sp->np_filler_depth, seed, csize.X, csize.Z);
	MapgenBasic::np_cave1 = sp->np_cave1;
	MapgenBasic::np_cave2 = sp->np_cave2;

	INodeDefManager *ndef = emerge->ndef;
	c_sand = ndef->getId("mapgen_sand");
}

MapgenCanyon::~MapgenCanyon()
{
	for(ChunkStorage::iterator it = m_chunks.begin(); it != m_chunks.end(); ++it) {
		delete it->second;
	}
}

MapgenCanyonParams::MapgenCanyonParams()
	: MapgenSpecificParams()
{
	spflags = 0;

	// This is the fixed value for 65kx65k maps and CanyonChunkContentSize=64.
	// It could be reduce for debugging the lower-res versions of the map
	world_level = 10 * 2;

	rain_factor = 0.0001;
	river_dominance_expansion = 1.3;
	river_fitting = 5.0;
	river_nice_routing = 1.5;
	river_rounding = 3;
	river_valley_broadening = 50.0;
	river_width_to_depth_ratio = 12.0;

	ground_base_scale = 1.4;
	ground_noise_scale = 0.15;
	ground_noise_off = 0.25;

	cave_width       = 0.2;

	np_filler_depth = NoiseParams(0, 1.2, v3f(150, 150, 150), 261,   3, 0.7, 2.0);
	np_cave1        = NoiseParams(0, 12,  v3f(61,  61,  61),  52534, 3, 0.5, 2.0);
	np_cave2        = NoiseParams(0, 12,  v3f(67,  67,  67),  10325, 3, 0.5, 2.0);
}

void MapgenCanyonParams::readParams(const Settings* settings)
{
	settings->getFlagStrNoEx("mgcanyon_spflags", spflags, flagdesc_mapgen_canyon);

	settings->getFloatNoEx("rain_factor", rain_factor);
	settings->getFloatNoEx("river_valley_broadening", river_valley_broadening);
	settings->getFloatNoEx("river_width_to_depth_ratio", river_width_to_depth_ratio);

	settings->getFloatNoEx("ground_base_scale", ground_base_scale);
	settings->getFloatNoEx("ground_noise_scale", ground_noise_scale);
	settings->getFloatNoEx("ground_noise_off", ground_noise_off);

	settings->getNoiseParams("mgcanyon_np_filler_depth", np_filler_depth);
	settings->getNoiseParams("mgcanyon_np_cave1",        np_cave1);
	settings->getNoiseParams("mgcanyon_np_cave2",        np_cave2);
}

void MapgenCanyonParams::writeParams(Settings* settings) const
{
	settings->setFlagStr("mgflat_spflags", spflags, flagdesc_mapgen_canyon, U32_MAX);

	settings->setFloat("rain_factor", rain_factor);
	settings->setFloat("river_valley_broadening", river_valley_broadening);
	settings->setFloat("river_width_to_depth_ratio", river_width_to_depth_ratio);

	settings->setFloat("ground_base_scale", ground_base_scale);
	settings->setFloat("ground_noise_scale", ground_noise_scale);
	settings->setFloat("ground_noise_off", ground_noise_off);

	settings->setNoiseParams("mgcanyon_np_filler_depth", np_filler_depth);
	settings->setNoiseParams("mgcanyon_np_cave1",        np_cave1);
	settings->setNoiseParams("mgcanyon_np_cave2",        np_cave2);
}

void MapgenCanyon::makeChunk(BlockMakeData* data)
{
	// Pre-conditions
	assert(data->vmanip);
	assert(data->nodedef);
	assert(data->blockpos_requested.X >= data->blockpos_min.X &&
		data->blockpos_requested.Y >= data->blockpos_min.Y &&
		data->blockpos_requested.Z >= data->blockpos_min.Z);
	assert(data->blockpos_requested.X <= data->blockpos_max.X &&
		data->blockpos_requested.Y <= data->blockpos_max.Y &&
		data->blockpos_requested.Z <= data->blockpos_max.Z);
	
	this->generating = true;
	this->vm = data->vmanip;
	this->ndef = data->nodedef;
	
	const v3s16 blockpos_min = data->blockpos_min;
	const v3s16 blockpos_max = data->blockpos_max;
	node_min = blockpos_min * MAP_BLOCKSIZE;
	node_max = (blockpos_max + 1) * MAP_BLOCKSIZE - 1;
	full_node_min = (blockpos_min - 1) * MAP_BLOCKSIZE;
	full_node_max = (blockpos_max + 2) * MAP_BLOCKSIZE - 1;
	
	blockseed = getBlockSeed2(full_node_min, seed);
	
	m_bgen->calcBiomeNoise(node_min);
	
	const v3s16 em = vm->m_area.getExtent();
	const MapNode n_air(CONTENT_AIR);
	const MapNode n_river_water(c_river_water_source);
	const MapNode n_sand(c_sand);
	const MapNode n_stone(c_stone);
	const MapNode n_water(c_water_source);
	const int row_size = node_max.X - node_min.X + 1;
	s16 stone_surface_max_y = -MAX_MAP_GENERATION_LIMIT;
	for(s16 z_chunk = node_min.Z; z_chunk <= node_max.Z; z_chunk += CanyonChunkContentSize)
	for(s16 x_chunk = node_min.X; x_chunk <= node_max.X; x_chunk += CanyonChunkContentSize) {
		const CanyonPatch *patch = getPatch(x_chunk, z_chunk, world_level);
		if(!patch)
			continue;

		const int x_min = std::max(0, node_min.X - patch->pos.X);
		const int z_min = std::max(0, node_min.Z - patch->pos.Y);
		const int x_max = std::min(node_max.X - patch->pos.X, CanyonChunkContentSize - 1);
		const int z_max = std::min(node_max.Z - patch->pos.Y, CanyonChunkContentSize - 1);
		for(s16 z = z_min; z <= z_max; z++)
		for(s16 x = x_min; x <= x_max; x++) {
			float l_ground = 0, l_river_bed = 0, l_river = 0;
			const int patchIndex = (x + CanyonChunkBorderSize) +
					(z + CanyonChunkBorderSize) * CanyonChunkBufferSize;
			const float xn = -patch->data[patchIndex].flow_x;
			const float xp = patch->data[patchIndex + 1].flow_x;
			const float yn = -patch->data[patchIndex].flow_y;
			const float yp = patch->data[patchIndex + CanyonChunkBufferSize].flow_y;
			l_ground = patch->data[patchIndex].ground_level;
			const float flow = std::max(0.0f, xn) + std::max(0.0f, yn) +
				std::max(0.0f, xp) + std::max(0.0f, yp);
			const float river_dominance = patch->data[patchIndex].river_dominance;
			if(flow > 1) {
				l_river_bed = 1;
				l_river = patch->data[patchIndex].fluid_level;
				l_ground = l_river - sqrt(flow) - l_river_bed;
			} else if(river_dominance > 0.999) {
				// Fake a river-bed beneath the river itself
				// TODO Limit river_dominance in vertical direction
				l_river_bed = 1;
			}
			l_river_bed += l_ground;

			if (l_ground > stone_surface_max_y)
				stone_surface_max_y = ceil(l_ground);

			u32 index_data = vm->m_area.index(patch->pos.X + x, node_min.Y - 1, patch->pos.Y + z);
			for(s16 y = node_min.Y - 1; y <= node_max.Y + 1; ++y) {
				if (vm->m_data[index_data].getContent() == CONTENT_IGNORE) {
					if(y < l_ground)
						vm->m_data[index_data] = n_stone;
					else if(y < l_river_bed)
						vm->m_data[index_data] = n_sand;
					else if(y < l_river)
						vm->m_data[index_data] = n_river_water;
					else
						vm->m_data[index_data] = n_air;
				}
				vm->m_area.add_y(em, index_data, 1);
			}
			
			u32 noise_index = (patch->pos.Y + z - node_min.Z) * row_size + patch->pos.X + x - node_min.X;
			float heat = m_bgen->heatmap[noise_index];
			heat = 50 + (1 - river_dominance) * (heat - 50);
			m_bgen->heatmap[noise_index] = heat;
			float humidity = m_bgen->humidmap[noise_index];
			humidity = 100 + (1 - river_dominance) * (humidity - 100);
			m_bgen->humidmap[noise_index] = humidity;
		}
	}
	
	updateHeightmap(node_min, node_max);
	
	MgStoneType stone_type = generateBiomes();
	
	// Cave creation.
	if (flags & MG_CAVES)
		generateCaves(stone_surface_max_y, -33);

	// Dungeon creation
	if ((flags & MG_DUNGEONS) && node_max.Y < 50)
		generateDungeons(stone_surface_max_y, stone_type);

	// Generate the registered decorations
	if (flags & MG_DECORATIONS)
		m_emerge->decomgr->placeAllDecos(this, blockseed, node_min, node_max);

	// Generate the registered decorations
	if (flags & MG_DECORATIONS)
		m_emerge->decomgr->placeAllDecos(this, blockseed, node_min, node_max);
	
	// Generate the registered ores
	m_emerge->oremgr->placeAllOres(this, blockseed, node_min, node_max);

	// Sprinkle some dust on top after everything else was generated
	dustTopNodes();

	// Add top and bottom side of water to transforming_liquid queue
	updateLiquid(&data->transforming_liquid, full_node_min, full_node_max);
	
	if (flags & MG_LIGHT)
		calcLighting(node_min - v3s16(0, 1, 0), node_max + v3s16(0, 1, 0),
					 full_node_min, full_node_max);
	
	this->generating = false;
}

int MapgenCanyon::getSpawnLevelAtPoint(v2s16 p)
{
	const CanyonPatch *patch = getPatch(p.X, p.Y, world_level);
	if(!patch) {
		return MAX_MAP_GENERATION_LIMIT;
	}

	// Decide spawn-level from patch
	const int patchIndex = (p.X - patch->pos.X + CanyonChunkBorderSize) +
			(p.Y - patch->pos.Y + CanyonChunkBorderSize) * CanyonChunkBufferSize;
	return patch->data[patchIndex].ground_level;
}

CanyonPatch* MapgenCanyon::getPatch(int x, int y, int level)
{
	const v2s32 map_pos = v2s32(x, y) + ((CanyonChunkContentSize / 2) << (level / 2));
	const int max_size = CanyonChunkContentSize << (level / 2);
	if(map_pos.X < 0 || map_pos.Y < 0 || map_pos.X > max_size || map_pos.Y > max_size)
		return NULL;
	CanyonPatch *result = getChunk(map_pos / CanyonChunkContentSize, level);
	return result;
}

CanyonPatch* MapgenCanyon::getChunk(v2s32 pos, int level)
{
	ChunkStorage::const_iterator it = m_chunks.find(v3s32(pos.X, pos.Y, level));
	if(it != m_chunks.end())
		return it->second;
	
	CanyonPatch * current = NULL;
	if(level > 0) {
		const CanyonPatch *const parent = getChunk(v2s32(pos.Y, pos.X / 2), level - 1);
		current = divide(parent, pos.X & 1);
	} else {
		current = createRoot();
	}
	m_chunks.insert(std::make_pair(v3s32(pos.X, pos.Y, level), current));
	return current;
}

CanyonPatch* MapgenCanyon::createRoot()
{
	CanyonPatch *root = new CanyonPatch();
	root->pos = v2s16(-CanyonChunkContentSize / 2, -CanyonChunkContentSize / 2);
	root->level = 0;
	for(int y = 0; y < CanyonChunkBufferSize; ++y) {
		const int p_y = y - CanyonChunkBorderSize;
		const float flow_y = (p_y - CanyonChunkContentSize / 2) * rain_factor;
		const int parent = p_y < CanyonChunkContentSize / 2 ? 1 : 3;
		int index = CanyonChunkBufferSize * y;
		for(int x = 0; x < CanyonChunkBufferSize; ++x) {
			CanyonNode &node = root->data[index];
			node.fluid_level = 0;
			node.ground_level = 0;
			node.flow_x = 0.0f;
			node.flow_y = flow_y;
			node.parent = parent;
			node.join_x = false;
			node.join_y = flow_y != 0;
			node.river_dominance = 0;
			++index;
		}
	}
	return root;
}

CanyonPatch *MapgenCanyon::divide(const CanyonPatch* parent, bool right_side)
{
	CanyonPatch *patch = new CanyonPatch();
	patch->pos = v2s16(parent->pos.Y * 2 + (right_side ? CanyonChunkContentSize : 0), parent->pos.X);
	patch->level = parent->level + 1;
	const int level = patch->level;
	const int scale_exp = (this->world_level - 1) / 2 - (level - 1) / 2;
	const bool isSecondSplit = (level - 1) & 0x1;
	const float broadening_limit = 1. / river_width_to_depth_ratio;
	// Increases on the last levels to form 1x1 nodes instead of e.g. 2x0.5
	// rivers which would result in solid ground interrupting the streams
	const float corrected_broadening_limit = std::max(double(broadening_limit),
			pow(0.25f, scale_exp));
	
	// Decide which child becomes the common exit
	// (false for negative side, true for positive)
	float splits[CanyonChunkBufferSize][CanyonChunkBufferSize / 2];
	for(int y = 0; y < CanyonChunkBufferSize; ++y) {
		int p_y = patch->pos.Y - CanyonChunkBorderSize + y;
		CanyonBaseCursor base(*parent, 0, y, right_side);
		for(int x = 0; x < CanyonChunkBufferSize / 2; ++x) {
			int p_x = patch->pos.X - CanyonChunkBorderSize + x * 2;
			const int parent = base.parent();
			if(parent & 1) {
				float prob_left = 1.0f, prob_right = 1.0f;
				
				if(base.flow() < corrected_broadening_limit) {
					// Cell consists of mostly solid land blocks
					
					const float x_neg = -base.flow_nx();
					const float x_pos = -base.flow_px();
					prob_left *= x_neg < 0.0 ? 1.0 :
								x_neg > 0.5 ? 0.0 :
								1 - 2 * x_neg;
					prob_right *= x_pos < 0.0 ? 1.0 :
								x_pos > 0.5 ? 0.0 :
								1 - 2 * x_pos;
								
					// Apply river_fitting
					if(base.at(-1, 0).ground_level() > base.at(1, 0).ground_level()) {
						prob_right *= river_fitting;
					} else {
						prob_left *= river_fitting;
					}
					
					// Apply river_route_nice to avoid cell borders
					prob_left *= pow(river_nice_routing, log(((p_x / 2 + 1) ^ (p_x / 2)) & 0xF));
					prob_right *= pow(river_nice_routing, log(((p_x / 2 - 1) ^ (p_x / 2)) & 0xF));
					
					// Apply some rounding of river corners
					int out_parent = -1;
					if(base.parent() == 1 && y > 0) {
						out_parent = base.at(0, -1).parent();
					} else if(base.parent() == 3 && y < CanyonChunkBufferSize - 1) {
						out_parent = base.at(0, 1).parent();
					}
					if(out_parent == 0) {
						prob_left *= river_rounding;
					} else if(out_parent == 2) {
						prob_right *= river_rounding;
					}
					
					const float noise = noise3d(p_x, p_y, level, seed) * 0.5 + 0.5;
					splits[y][x] = (noise * (prob_left + prob_right) - prob_left > 0) ? 1 : 0;
				} else {
					// Cell consists of mostly river blocks
					
					float self_level = base.fluid_level();
					float neg_level = base.at(-1,0).fluid_level();
					float pos_level = base.at(1,0).fluid_level();
					if(neg_level < pos_level && neg_level < self_level && !base.join_nx()) {
						// Avoid a cliff on the -X side
						splits[y][x] = 1.0;
					} else if(pos_level < self_level && !base.join_px()) {
						// Avoid a cliff on the +X side
						splits[y][x] = 0.0;
					} else {
						// Tries to concentrate the main flow of the river in its center
						float flow_y = base.flow_y();
						float prob_left = std::max(0.0f, (flow_y + base.at(-1, 0).flow_y()) / flow_y - 0.2f);
						float prob_right = std::max(0.0f, (flow_y + base.at(1, 0).flow_y()) / flow_y - 0.2f);
						float left = std::max(0.0, prob_left * (noise3d(p_x, p_y, level, seed + 47) * 0.5 + 0.5) - 0.2);
						float right = std::max(0.0, prob_right * (noise3d(p_x, p_y, level, seed+ 49) * 0.5 + 0.5) - 0.2);
						if(left + right == 0) {
							splits[y][x] = noise3d(p_x, p_y, level, seed + 42) * 0.5 + 0.5;
						} else {
							splits[y][x] = right / (left + right);
						}
					
						// Reduce irregular splitting where possible
						float smooth_factor = std::max(0.0, base.flow() - 1.0);
						splits[y][x] = (1 * splits[y][x] + smooth_factor * 0.5) / (1 + smooth_factor);
					}
				}
			} else {
				// Decide split direction based on main flow exit
				splits[y][x] = (parent & 2) ? 1 : 0;
			}
			base.move(1, 0);
		}
	}
	
	// Initialize child nodes' flow at borders of old node
	for(int y = 1; y < CanyonChunkBufferSize; ++y) {
		CanyonBaseCursor base(*parent, 0, y, right_side);
		CanyonEditCursor edit(*patch, 0, y);
		for(int x = 0; x < CanyonChunkBufferSize / 2; ++x) {
			const float split_prev = splits[y-1][x];
			const float split = splits[y][x];
			
			// Split the ny border
			edit.set_left_flow_nx(base.flow_nx());
			edit.set_left_join_nx(base.join_nx());
			edit.set_right_join_nx(true);
			if(base.at(0, -1).parent() == 3) {
				edit.set_left_flow_ny((1 - split_prev) * base.flow_ny());
				edit.set_right_flow_ny(split_prev * base.flow_ny());
				edit.set_left_join_ny(split_prev < 0.99);
				edit.set_right_join_ny(split_prev > 0.01);
			} else if(base.parent() == 1) {
				edit.set_left_flow_ny((1 - split) * base.flow_ny());
				edit.set_right_flow_ny(split * base.flow_ny());
				edit.set_left_join_ny(split < 0.99);
				edit.set_right_join_ny(split > 0.01);
			} else {
				edit.set_left_flow_ny(0.5 * base.flow_ny());
				edit.set_right_flow_ny(0.5 * base.flow_ny());
				edit.set_left_join_ny(base.flow_ny() > 0.01);
				edit.set_right_join_ny(base.flow_ny() > 0.01);
			}
			
			base.move(1, 0);
			edit.move(1, 0);
		}
	}
	
	// Initialize child nodes' flow at internal border of old node
	for(int y = 1; y < CanyonChunkBufferSize - 1; ++y) {
		CanyonEditCursor edit(*patch, 0, y);
		for(int x = 0; x < CanyonChunkBufferSize / 2; ++x) {
			float left_flow = edit.left_flow_nx() + edit.left_flow_ny() +
					edit.left_flow_py();
			edit.set_right_flow_nx(left_flow + (isSecondSplit ? 0.25 : 0.5) * rain_factor);

			edit.move(1, 0);
		}
	}
	
	// Initialize fluid_level and parent-direction
	for(int y = 1; y < CanyonChunkBufferSize - 1; ++y) {
		int p_y = patch->pos.Y - CanyonChunkBorderSize + y;
		CanyonBaseCursor base(*parent, 1, y, right_side);
		CanyonEditCursor edit(*patch, 1, y);
		for(int x = 1; x < CanyonChunkBufferSize / 2 - 1; ++x) {
			int p_x = patch->pos.X - CanyonChunkBorderSize + x * 2;
			const bool split = edit.right_flow_nx() > 0;
			
			float min_level = base.fluid_level();
			float max_level = base.ground_level() * 1;
			
 			const float level_noise = noise3d(p_x, p_y, level, seed + 78) * 0.5 + 0.5;
			const float level_noise_adj = pow(level_noise , 0.3);
			if(split) {
				if(edit.left_flow_nx() > 0)
					max_level = std::min(base.at(-1, 0).fluid_level(), max_level);
				if(edit.left_flow_ny() > 0)
					max_level = std::min(base.at(0, -1).fluid_level(), max_level);
				if(edit.left_flow_py() > 0)
					max_level = std::min(base.at(0, 1).fluid_level(), max_level);
				max_level = std::max(min_level, max_level);
				edit.set_left_fluid_level(min_level + (max_level - min_level) * level_noise_adj);
				edit.set_right_fluid_level(min_level);
			} else {
				if(edit.right_flow_px() > 0)
					max_level = std::min(base.at(1, 0).fluid_level(), max_level);
				if(edit.right_flow_ny() > 0)
					max_level = std::min(base.at(0, -1).fluid_level(), max_level);
				if(edit.right_flow_py() > 0)
					max_level = std::min(base.at(0, 1).fluid_level(), max_level);
				max_level = std::max(min_level, max_level);
				edit.set_left_fluid_level(min_level);
				edit.set_right_fluid_level(min_level + (max_level - min_level) * level_noise_adj);
			}
			edit.set_left_fluid_level(edit.left_fluid_level() * ground_base_scale);
			edit.set_right_fluid_level(edit.right_fluid_level() * ground_base_scale);
			
			if(split) {
				// TODO Allow different parent
				edit.set_left_parent(2);
				edit.set_right_parent(base.parent());
			} else {
				edit.set_left_parent(base.parent());
				// TODO Allow different parent
				edit.set_right_parent(0);
			}
			
			base.move(1, 0);
			edit.move(1, 0);
		}
	}
	
	// Update ground heights
	for(int y = 0; y < CanyonChunkBufferSize; ++y) {
		int p_y = patch->pos.Y - CanyonChunkBorderSize + y;
		CanyonBaseCursor base(*parent, 1, y, right_side);
		CanyonEditCursor edit(*patch, 1, y);
		for(int x = 1; x < CanyonChunkBufferSize / 2 - 1; ++x) {
			int p_x = patch->pos.X - CanyonChunkBorderSize + x * 2;
			const float near_flow = std::max(std::max(std::max(base.at(-1,0).flow(),
					base.at(1,0).flow()), base.at(0,-1).flow()), base.at(0,1).flow());
			const float flow = (2 * base.flow() + near_flow) / 3;
			
			// TODO Make this dependent on broadening_limit
			const float base_level = base.ground_level() +
					(base.fluid_level() - base.ground_level()) *
					std::min(flow * river_valley_broadening, 1.f);
			
			const float left_base = (base_level * 3 + base.at(-1, 0).ground_level() * 1) / 4;
			const float left_level = groundNoise(left_base, p_x, p_y, level);
			edit.set_left_ground_level(std::max(edit.left_fluid_level(), left_level));
// 			edit.set_left_ground_level(left_level);
			const float right_base = (base_level * 3 + base.at(1, 0).ground_level() * 1) / 4;
			const float right_level = groundNoise(right_base, p_x + 1, p_y, level);
			edit.set_right_ground_level(std::max(edit.right_fluid_level(), right_level));
// 			edit.set_right_ground_level(right_level);
			
			const float left_river = (base.river_dominance() * 3 + base.at(-1, 0).river_dominance() * 1) / 4;
			edit.set_left_river_dominance(std::min(std::max(flow, left_river * river_dominance_expansion), 1.0f));
			const float right_river = (base.river_dominance() * 3 + base.at(+1, 0).river_dominance() * 1) / 4;
			edit.set_right_river_dominance(std::min(std::max(flow, right_river * river_dominance_expansion), 1.0f));
			
			base.move(1, 0);
			edit.move(1, 0);
		}
	}
	
	// Randomly reduce fluid-level where this is possible
	for(int y = 0; y < CanyonChunkBufferSize; ++y) {
		int p_y = patch->pos.Y - CanyonChunkBorderSize + y;
		CanyonEditCursor edit(*patch, 1, y);
		for(int x = 1; x < CanyonChunkBufferSize / 2 - 1; ++x) {
			int p_x = patch->pos.X - CanyonChunkBorderSize + x * 2;
			
			float left_lowest = 0;
			switch(edit.left_parent()) {
				case 0:
					left_lowest = edit.at(-1, 0).right_fluid_level();
					break;
				case 1:
					left_lowest = edit.at(0, -1).left_fluid_level();
					break;
				case 2:
					left_lowest = edit.right_fluid_level();
					break;
				case 3:
					left_lowest = edit.at(0, 1).left_fluid_level();
					break;
			}
			const float left_noise = pow(noise3d(p_x, p_y, level, seed + 101) * 0.5 + 0.5, 0.2);
			edit.set_left_fluid_level(left_noise * (edit.left_fluid_level() - left_lowest) + left_lowest);
			
			float right_lowest = 0;
			switch(edit.right_parent()) {
				case 0:
					right_lowest = edit.left_fluid_level();
					break;
				case 1:
					right_lowest = edit.at(0, -1).right_fluid_level();
					break;
				case 2:
					right_lowest = edit.at(1, 0).left_fluid_level();
					break;
				case 3:
					right_lowest = edit.at(0, 1).right_fluid_level();
					break;
			}
			const float right_noise = pow(noise3d(p_x, p_y, level, seed + 102) * 0.5 + 0.5, 0.2);
			edit.set_right_fluid_level(right_noise * (edit.right_fluid_level() - right_lowest) + right_lowest);

			edit.move(1, 0);
		}
	}
	
	// Scale flow up with the same factor as the landscape
	if(isSecondSplit) {
		for(int y = 0; y < CanyonChunkBufferSize; ++y) {
			CanyonEditCursor edit(*patch, 0, y);
			for(int x = 0; x < CanyonChunkBufferSize / 2; ++x) {
				edit.set_left_flow_nx(edit.left_flow_nx() * 4);
				edit.set_left_flow_ny(edit.left_flow_ny() * 4);
				edit.set_right_flow_nx(edit.right_flow_nx() * 4);
				edit.set_right_flow_ny(edit.right_flow_ny() * 4);
				edit.move(1, 0);
			}
		}
	}
	
	return patch;
}

float MapgenCanyon::groundNoise(float base, int x, int y, int level)
{
	const float noise = noise3d(x, y, level, seed) + ground_noise_off;
	return base * ground_base_scale + ground_noise_scale * noise;
}
