#ifndef PAIRHASH_H
#define PAIRHASH_H

#include <utility>
#include <cstdint>

struct PairHash
{

	static uint32_t hash_combine(uint32_t lhs, uint32_t rhs) {
		lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
		return lhs;
	}

	static uint64_t hash_combine(uint64_t lhs, uint64_t rhs) {
		lhs ^= rhs + 0x9e3779b97f4a7c15 + (lhs << 6) + (lhs >> 2);
		return lhs;
	}

	template <class T1, class T2>
	std::size_t operator() (std::pair<T1, T2> const &pair) const
	{
		std::size_t h1 = std::hash<T1>()(pair.first);
		std::size_t h2 = std::hash<T2>()(pair.second);
		return hash_combine(h1, h2);

//		size_t hash_combine( size_t lhs, size_t rhs ) {
//			lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
//			return lhs;
//		}

	}
};

#endif // PAIRHASH_H
