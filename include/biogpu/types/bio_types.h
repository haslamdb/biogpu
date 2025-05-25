#ifndef BIOGPU_BIO_TYPES_H
#define BIOGPU_BIO_TYPES_H

#include <cstdint>
#include <string>

namespace biogpu {

// Forward declarations
struct DNASequence;
struct KmerIndex;
struct AlignmentResult;
struct ResistanceMutation;

// DNA encoding
enum class DNABase : uint8_t {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4  // Unknown
};

// Quality score type
using QualityScore = uint8_t;

// Sequence ID type
using SequenceID = uint32_t;

} // namespace biogpu

#endif // BIOGPU_BIO_TYPES_H
