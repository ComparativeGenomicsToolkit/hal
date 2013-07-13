#ifndef _COPY_GENOMES_H
#define _COPY_GENOMES_H
namespace hal {
class Genome;
}
// copy genome dimensions (sequence size, top segments size, bot segments size)
void copyAllDimensions(hal::AlignmentConstPtr inAlignment,
                       const hal::Genome *inGenome, hal::Genome *outGenome);
// copy only bottom segment dimensions from one genome to another
void copyBotDimensions(const hal::Genome *inGenome, hal::Genome *outGenome);
// copy DNA, top segments, metadata from one genome to another
void copyGenomeWithoutBotSegments(const hal::Genome *inGenome,
                                  hal::Genome *outGenome);
// copy bottom segments from one genome to another, keeping segments
// pointing to the proper children
void copyBotSegments(const hal::Genome *inGenome, hal::Genome *outGenome);
// fix the bottom/top parse info
void fixParseInfo(hal::Genome *genome);
// copy only top segment dimensions from one genome to another
void copyTopDimensions(const hal::Genome *inGenome, hal::Genome *outGenome);
#endif // _COPY_GENOMES_H
