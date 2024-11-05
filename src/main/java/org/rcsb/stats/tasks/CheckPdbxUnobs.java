package org.rcsb.stats.tasks;

import org.rcsb.cif.CifIO;
import org.rcsb.cif.ParsingException;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.*;
import org.rcsb.stats.Constants;
import org.rcsb.stats.Helpers;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Obtain the number of non-hydrogen atoms as described by the `atom_site` content of an mmCIF file.
 */
public class CheckPdbxUnobs {
    static final Logger logger = LoggerFactory.getLogger(CheckPdbxUnobs.class);
    final AtomicInteger counter = new AtomicInteger();

    public static void main(String[] args) throws IOException {
        CheckPdbxUnobs me = new CheckPdbxUnobs();
        me.countHeavyAtoms(me.fetchStructureData("1IBA"));
        new CheckPdbxUnobs().computeStats();
    }

    void computeStats() throws IOException {
        // obtain set of all known identifiers (1ABC, 1ABD, ...)
        Set<String> identifiers = Helpers.getAllIdentifiers(Constants.ResultsContentType.EXPERIMENTAL);

        // traverse all identifiers
        long heavyAtomCount = identifiers.parallelStream().peek(this::logProgress)
                // load the {@link CifFile} for each identifier
                .map(this::fetchStructureData)
                // process each structure: count the number of non-hydrogen atoms
                .mapToLong(this::countHeavyAtoms)
                // aggregate as sum
                .sum();

        logger.info("There are {} heavy (non-hydrogen) atoms in {} PDB structures", Helpers.formatNumber(heavyAtomCount), Helpers.formatNumber(counter.get()));
    }

    /**
     * Obtain structure for a single identifier from models.rcsb.org.
     * @param identifier datum to load
     * @return structure data in the form of a {@link CifFile}
     */
    CifFile fetchStructureData(String identifier) {
        try {
            URL url = new URL(String.format(Constants.BCIF_SOURCE, identifier));
            // other CifIO methods allow reading from Paths or byte streams -- methods for writing can be found there too
            return CifIO.readFromURL(url);
        } catch (IOException e) {
            logger.warn("Failed to pull structure data for {}", identifier);
            throw new UncheckedIOException(e);
        } catch (ParsingException e) {
            logger.warn("Failed to parse structure data for {}", identifier);
            throw e;
        }
    }

    /**
     * Process a CIF file (obtained from @link{CifIO#readFromUrl()}).
     * @param cifFile source data
     * @return the count of non-hydrogen atoms
     */
    long countHeavyAtoms(CifFile cifFile) {
        Map<String, Integer> entityIdToLength = new HashMap<>();
        Map<String, String> asymIdsToEntId  = new HashMap<>();
        MmCifBlock block = cifFile
                // optional: apply mmCIF schema to get schema definitions and types
                .as(StandardSchemata.MMCIF)
                // CIF files may have multiple blocks of data, the PDB archive only makes use of the 1st
                .getFirstBlock();

        String entryId = block.getStruct().getEntryId().values().findFirst().orElse("");

        PdbxNmrEnsemble ensemble = block.getPdbxNmrEnsemble();
        int numModels = ensemble.getConformersSubmittedTotalNumber().values().findAny().orElse(1);

        EntityPolySeq entityPoly = block.getEntityPolySeq();
        for (int rowIndex = 0; rowIndex < entityPoly.getRowCount(); rowIndex++) {
            String entId = entityPoly.getEntityId().get(rowIndex);
            entityIdToLength.merge(entId, 1, Integer::sum);
        }

        AtomSite atomSite = block.getAtomSite();
        Map<String, TreeSet<Integer>> asymIdToResNums = new HashMap<>();
        double totalResOcc = 0.0;
        int prevSeqId = -1;
        String prevAsymId = null;
        for (int rowIndex = 0; rowIndex < atomSite.getRowCount(); rowIndex++) {
            String entId = atomSite.getLabelEntityId().get(rowIndex);
            String asymId = atomSite.getLabelAsymId().get(rowIndex);
            asymIdsToEntId.putIfAbsent(asymId, entId);
            if (!entityIdToLength.containsKey(entId)) continue; // not a polymer
            int seqId = atomSite.getLabelSeqId().get(rowIndex);
            double occ = atomSite.getOccupancy().get(rowIndex);

            if (prevSeqId>0 && seqId!=prevSeqId && totalResOcc > 0.00001) {
                asymIdToResNums.computeIfAbsent(prevAsymId, k -> new TreeSet<>()).add(prevSeqId);
            }

            if (prevSeqId>0 && seqId!=prevSeqId) {
                totalResOcc = 0.0;
            }
            totalResOcc += occ;

            prevSeqId = seqId;
            prevAsymId = asymId;
        }
        if (totalResOcc > 0.00001) {
            asymIdToResNums.computeIfAbsent(prevAsymId, k -> new TreeSet<>()).add(prevSeqId);
        }

        Map<String, List<String>> entIdToAsymIds = new HashMap<>();
        asymIdsToEntId.entrySet().stream().forEach(e -> entIdToAsymIds.computeIfAbsent(e.getValue(), k -> new ArrayList<>()).add(e.getKey()));

        Map<String, Integer> asymIdToUnmodeledCount = new HashMap<>();
        for (String entId : entityIdToLength.keySet()) {
            for (String asymId : entIdToAsymIds.get(entId)) {
                TreeSet<Integer> allResNums = asymIdToResNums.get(asymId);
                if (allResNums == null) {
                    logger.info("No residue numbers info present for entry {}, asym {}", entryId, asymId);
                    continue;
                }
                for (int i = 1; i <= entityIdToLength.get(entId); i++) {
                    if (!allResNums.contains(i)) {
                        asymIdToUnmodeledCount.merge(asymId, 1, Integer::sum);
                    }
                }
            }
        }

        PdbxUnobsOrZeroOccResidues pdbxUnobsRes = block.getPdbxUnobsOrZeroOccResidues();
        Map<String, Integer> asymIdToPdbUnobsCount = new HashMap<>();
        for (int rowIndex = 0; rowIndex < pdbxUnobsRes.getRowCount(); rowIndex++) {
            String asymId = pdbxUnobsRes.getLabelAsymId().get(rowIndex);
            asymIdToPdbUnobsCount.merge(asymId, 1, Integer::sum);
        }

        for (String entId : entityIdToLength.keySet()) {
            for (String asymId : entIdToAsymIds.get(entId)) {
                int unmodeledCount = Optional.ofNullable(asymIdToUnmodeledCount.get(asymId)).orElse(0) * numModels;
                int unobsCount = Optional.ofNullable(asymIdToPdbUnobsCount.get(asymId)).orElse(0);
                if (unmodeledCount != unobsCount) {
                    logger.info("Different count for entry {}, asym {}, unmodeled {}, unobs {}", entryId, asymId, unmodeledCount, unobsCount);
                }
            }
        }

//        long countUnobsRes = block.getPdbxUnobsOrZeroOccResidues()
//                .getLabelSeqId()
//                .values().count();
//
//        if (countUnobsRes > 0) {
//            logger.info("Entry {}. Unobserved residues from pdbx_unobs_or_zero_occ_residues: {}", block.getStruct().getEntryId().values().findFirst().orElse(""), countUnobsRes);
//        }

//        long countUnobsAtom = block.getPdbxUnobsOrZeroOccAtoms()
//                .getLabelSeqId()
//                .values().count();
//        if (countUnobsAtom > 0) {
//            logger.info("Entry {}. Unobserved atoms from pdbx_unobs_or_zero_occ_atom: {}", block.getStruct().getEntryId().values().findFirst().orElse(""), countUnobsAtom);
//        }

        return 0;
    }

    final Set<String> HYDROGEN_ATOMS = Set.of("H", "D", "T");

    /**
     * Filters for non-hydrogen atoms based on their `atom_site.type_symbol`.
     * @param typeSymbol element of this atom
     * @return false if this is hydrogen
     */
    boolean isHeavyAtom(String typeSymbol) {
        return !HYDROGEN_ATOMS.contains(typeSymbol);
    }

    /**
     * Indicate that 10,000 entries have been processed.
     */
    void logProgress(Object ignored) {
        if (counter.incrementAndGet() % 10000 == 0) {
            logger.info("Processed {} entries", Helpers.formatNumber(counter.get()));
        }
    }
}
