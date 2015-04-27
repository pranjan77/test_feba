
CREATE TABLE Organism(
   orgId            TEXT     NOT NULL,
   division         TEXT     NOT NULL,
   genus            TEXT     NOT NULL,
   species          TEXT     NOT NULL,
   strain           TEXT     NOT NULL,
   taxonomyId       INT,     /* NCBI taxonomyId */
   PRIMARY KEY (orgId)
);

CREATE TABLE Gene(
   orgId            TEXT     NOT NULL,
   locusId          TEXT     NOT NULL, /* nickname x locusId is unique; locusId may not be */
   sysName          TEXT,    /* a locus tag like SO_1446 or b2338, sometimes identical to locusId */
   scaffoldId       TEXT     NOT NULL,
   begin            INT      NOT NULL,
   end              INT      NOT NULL,
   /* Type is: 1 for protein-coding, 2 for rRNA, 5 for tRNA, 6 for ncRNA, 7 for pseudogene,
      9 for CRISPR repeat, 10 for CRISPR spacer, 11 for antisense RNA, 99 for other (possibly a pseudogene)
   */
   type             INT      NOT NULL,
   strand           TEXT     NOT NULL,
   gene             TEXT,    /* a gene name like recA */
   desc             TEXT,
   GC               REAL,    /* %GC of the gene's sequence */
   nTA              INT,     /* number of TA dinucleotides in the gene's sequence */
   PRIMARY KEY (orgId, locusId)
);

/* Only experiments that succeeded should be loaded into the database.
   Quality metrics are documented in lib/FEBA_template.html
 */
CREATE TABLE Experiment(
   orgId         TEXT       NOT NULL,
   expName       TEXT       NOT NULL,
   expDesc       TEXT       NOT NULL, /* the short form */
   timeZeroSet   TEXT       NOT NULL, /* which set of Time0 samples were used as a control */
   num           INT        NOT NULL, /* a secondary identifier, unique within each organism */
   nMapped       INT        NOT NULL,
   nPastEnd      INT        NOT NULL,
   nGenic        INT        NOT NULL,
   nUsed         INT        NOT NULL,
   gMed          INT        NOT NULL,
   gMedt0        INT        NOT NULL,
   gMean         REAL       NOT NULL,
   cor12         REAL       NOT NULL,
   mad12         REAL       NOT NULL,
   mad12c        REAL       NOT NULL,
   mad12c_t0     REAL       NOT NULL,
   opcor         REAL       NOT NULL,
   adjcor        REAL       NOT NULL,
   gccor         REAL       NOT NULL,
   maxFit        REAL       NOT NULL,
   PRIMARY KEY (orgId, expName)
);

CREATE TABLE GeneFitness(
   orgId            TEXT     NOT NULL,
   locusId          TEXT     NOT NULL,
   expName          TEXT     NOT NULL,
   fit              REAL     NOT NULL,
   t                REAL     NOT NULL,
   PRIMARY KEY (orgId,locusId,expName)
);

/* Most cofit genes for each gene.
   Genes in organisms that have relatively few experiments will not be included.
*/
CREATE TABLE Cofit(
	orgId TEXT NOT NULL,
        locusId TEXT NOT NULL,
	hitId TEXT NOT NULL,
        rank INT NOT NULL,
        cofit REAL NOT NULL,
        PRIMARY KEY (orgId,locusId,hitId)
);
