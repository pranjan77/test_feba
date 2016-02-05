
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

CREATE TABLE Ortholog(
   orgId1           TEXT     NOT NULL,
   locusId1         TEXT     NOT NULL,
   orgId2           TEXT     NOT NULL,
   locusId2         TEXT     NOT NULL,
   ratio            REAL     NOT NULL, /* BLAST_bit_score(alignment) / score(self_alignment) */
   PRIMARY KEY (orgId1,locusId1,orgId2,locusId2)
);

/* Only experiments that succeeded should be loaded into the database.
   Quality metrics are documented in lib/FEBA_template.html
 */
CREATE TABLE Experiment(
   orgId         TEXT       NOT NULL,
   expName       TEXT       NOT NULL, /* orgId x expName should be unique; expName may not be */
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
   expGroup      TEXT       NOT NULL,
   expDescLong   TEXT       NOT NULL,
   mutantLibrary TEXT       NOT NULL,
   person        TEXT       NOT NULL,
   dateStarted   TEXT       NOT NULL,
   setName       TEXT       NOT NULL,
   seqindex      TEXT       NOT NULL,
   media         TEXT       NOT NULL,
   /* These fields may be absent, but this should be represented as empty strings */
   temperature   TEXT       NOT NULL, /* should be in celsius */
   pH            TEXT       NOT NULL,
   vessel        TEXT       NOT NULL, /* Growth.Method in R tables */
   aerobic       TEXT       NOT NULL, /* Aerobic_v_Anaerobic */
   liquid        TEXT       NOT NULL, /* Liquid.v..solid */
   shaking       TEXT       NOT NULL,
   condition_1   TEXT       NOT NULL,
   units_1       TEXT       NOT NULL,
   concentration_1 TEXT     NOT NULL,
   condition_2   TEXT       NOT NULL,
   units_2       TEXT       NOT NULL,
   concentration_2 TEXT     NOT NULL,
   growthPlate TEXT NOT NULL,
   growthWells TEXT NOT NULL,
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

/* There is also, for each organism, a table
   FitByExp_org with fields expName,locusId,fit,t
   It is stored by experiment so it is fast to look
   up all the data for an experiment this way.
*/

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

/* Specific phenotypes -- a very sparse subset of gene/experiment combinations */
CREATE TABLE SpecificPhenotype(
	orgId TEXT NOT NULL,
	expName TEXT NOT NULL,
	locusId TEXT NOT NULL,
	PRIMARY KEY (orgId,expName,locusId)
);

CREATE TABLE GeneDomain(
   domainDb TEXT NOT NULL,
   orgId TEXT NOT NULL,
   locusId TEXT NOT NULL,
   domainId TEXT NOT NULL,
   domainName TEXT NOT NULL,
   begin INT NOT NULL,
   end INT NOT NULL,
   score REAL NOT NULL,
   evalue REAL NOT NULL,
   type TEXT,
   geneSymbol TEXT,
   ec TEXT,
   definition TEXT,
   /* Combination of domain, begin, end should be unique for each locus */
   PRIMARY KEY (orgId,locusId,domainId,begin,end)
);
CREATE INDEX 'orgLocus' on GeneDomain ('orgId' ASC, 'locusId' ASC);
CREATE INDEX 'domainDbId' on GeneDomain ('domainDb' ASC, 'domainId' ASC);
CREATE INDEX 'domainDbName' on GeneDomain ('domainDb' ASC, 'domainName' ASC);
CREATE INDEX 'domain_ec' on GeneDomain ('ec' ASC, 'orgId' ASC);

/* For each kilobase on each scaffold, shows the seek position into db.StrainFitness.orgId
   This arrangement is used because the StrainFitness tables are so large.
   If there are no insertions within 1000*kb to 1000*kb+999, then the entry might be omitted.
 */
CREATE TABLE StrainDataSeek(
       orgId TEXT NOT NULL,
       scaffoldId TEXT NOT NULL,
       kb INT NOT NULL,
       seek INT NOT NULL,
       PRIMARY KEY (orgId,scaffoldId,kb)
);

CREATE TABLE Compounds(
	compound TEXT NOT NULL,
        MW REAL,                /* molecular weight in g/mol */
        CAS TEXT,               /* CAS number usable at commonchemistry.org */
        PRIMARY KEY (compound)
);

CREATE TABLE MediaComponents(
	media TEXT NOT NULL,
        compound TEXT NOT NULL,
        concentration REAL,
        units TEXT
);
CREATE INDEX 'MediaComponentsByMedia' on MediaComponents('media' ASC);

CREATE TABLE LocusXref(
       orgId TEXT NOT NULL,
       locusId TEXT NOT NULL,
       xrefDb TEXT NOT NULL,
       xrefId TEXT NOT NULL);
CREATE INDEX 'LocusXrefByLocus' on LocusXref ('orgId', 'locusId', 'xrefDb');

/* For each protein in our genomes, its best hit in KEGG, if a good one exists */
CREATE TABLE BestHitKEGG(
       orgId TEXT NOT NULL,
       locusId TEXT NOT NULL,
       keggOrg TEXT NOT NULL, /* i.e., ajs from ajs:Ajs_1688 */
       keggId TEXT NOT NULL, /* i.e. Ajs_1688 */
       identity REAL NOT NULL, /* i.e., 99.0% amino acid identity */
       PRIMARY KEY (orgId,locusId)
);
CREATE INDEX 'KEGGToHits' on BestHitKEGG (keggOrg,keggId,orgId,locusId);

/* A KEGG gene may belong to more than one kgroup. Only information about hits is stored */
CREATE TABLE KEGGMember(
	keggOrg TEXT NOT NULL,
        keggId TEXT NOT NULL,
        kgroup TEXT NOT NULL, /* like K13522 */
        PRIMARY KEY (keggOrg,keggId,kgroup)
);
CREATE INDEX 'KEGGMemberByGroup' on KEGGMember (kgroup,keggOrg,keggId);

/* Description of each KEGG orthology group */
CREATE TABLE KgroupDesc(
       kgroup TEXT NOT NULL,
       desc TEXT NOT NULL,
       PRIMARY KEY (kgroup)
);

/* A kgroup may have more than one EC number. Only information relevant to hits is stored */
CREATE TABLE KgroupEC(
       kgroup TEXT NOT NULL,
       ecnum TEXT NOT NULL,
       PRIMARY KEY (kgroup,ecnum)
);
CREATE INDEX 'KgroupECByEcnum' on KgroupEC ('ecnum', 'kgroup');

/* An auxiliary table to speed up ortholog conditions page, with
   pre-selected ortholog groups of genes that have a specific phenotype
   in each condition and are BBHs of each other. Ideally, all genes in
   an OG would be BBHs of each other, but in practice, these are
   computed by clustering the BBHs.
*/
CREATE TABLE SpecOG(
       ogId INT NOT NULL,       /* arbitrary, unique for each ortholog group */
       expGroup TEXT NOT NULL,
       condition TEXT NOT NULL,
       orgId TEXT NOT NULL,
       locusId TEXT NOT NULL,
       minFit REAL NOT NULL,
       maxFit REAL NOT NULL,
       minT REAL NOT NULL,
       maxT REAL NOT NULL,
       PRIMARY KEY (expGroup, condition, orgId, locusId)
);
CREATE INDEX 'SpecOGByLocus' on SpecOG ('orgId','locusId');
CREATE INDEX 'SpecOGByOG' on SpecOG ('ogId');
