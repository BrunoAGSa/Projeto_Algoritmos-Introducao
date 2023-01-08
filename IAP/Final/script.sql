-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `mydb` DEFAULT CHARACTER SET utf8 ;
USE `mydb` ;

-- -----------------------------------------------------
-- Table `mydb`.`Locus`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Locus` (
  `ID` VARCHAR(50) NOT NULL,
  `Size` INT NOT NULL,
  `Molecular_type` VARCHAR(50) NOT NULL,
  `Genbank_division` VARCHAR(3) NOT NULL,
  `Modification_date` DATE NOT NULL,
  PRIMARY KEY (`ID`));


-- -----------------------------------------------------
-- Table `mydb`.`Sequences`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Sequences` (
  `ID` VARCHAR(50) NOT NULL,
  `Sequence` VARCHAR(2000) NOT NULL,
  `Length` INT NULL,
  `Count_A` INT NULL,
  `Count_C` INT NULL,
  `Count_T` INT NULL,
  `Count_G` INT NULL,
  PRIMARY KEY (`ID`));


-- -----------------------------------------------------
-- Table `mydb`.`References`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`References` (
  `ID` VARCHAR(50) NOT NULL,
  `Pubmed` VARCHAR(50) NOT NULL,
  `Title` VARCHAR(200) NULL,
  `Journal` VARCHAR(200) NULL,
  `Consortium` VARCHAR(200) NULL,
  `Remark` VARCHAR(200) NULL,
  PRIMARY KEY (`ID`, `Pubmed`));


-- -----------------------------------------------------
-- Table `mydb`.`Source`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Source` (
  `ID` VARCHAR(50) NOT NULL,
  `Location_span` VARCHAR(50) NOT NULL,
  PRIMARY KEY (`ID`));


-- -----------------------------------------------------
-- Table `mydb`.`CDS`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`CDS` (
  `ID_CDS` VARCHAR(50) NOT NULL,
  `Translacion_seq` VARCHAR(2000) NOT NULL,
  `Localization` VARCHAR(20) NOT NULL,
  `ID_protein` VARCHAR(50) NOT NULL,
  `Protein` VARCHAR(50) NOT NULL,
  PRIMARY KEY (`ID_CDS`));


-- -----------------------------------------------------
-- Table `mydb`.`Features`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Features` (
  `ID` VARCHAR(50) NOT NULL,
  `CDS` VARCHAR(50) NOT NULL,
  `Gene` INT NULL,
  `Regulatory` INT NULL,
  `exons` INT NULL,
  `poli_A-site` INT NULL,
  `misc_feature` INT NULL,
  `mRNA` INT NULL,
  PRIMARY KEY (`ID`),
  INDEX `features-CDS_idx` (`CDS` ASC) VISIBLE,
  CONSTRAINT `Features-Source`
    FOREIGN KEY (`ID`)
    REFERENCES `mydb`.`Source` (`ID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `features-CDS`
    FOREIGN KEY (`CDS`)
    REFERENCES `mydb`.`CDS` (`ID_CDS`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION);


-- -----------------------------------------------------
-- Table `mydb`.`Genbank`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Genbank` (
  `ID` VARCHAR(50) NOT NULL,
  `Organism` VARCHAR(200) NOT NULL,
  `Accession` VARCHAR(50) NOT NULL,
  `Definition` VARCHAR(200) NULL,
  `Keywords` VARCHAR(200) NULL,
  PRIMARY KEY (`ID`),
  CONSTRAINT `Genbank-locus`
    FOREIGN KEY (`ID`)
    REFERENCES `mydb`.`Locus` (`ID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Genbank-sequence`
    FOREIGN KEY (`ID`)
    REFERENCES `mydb`.`Sequences` (`ID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Genbank-references`
    FOREIGN KEY (`ID`)
    REFERENCES `mydb`.`References` (`ID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Genbank-features`
    FOREIGN KEY (`ID`)
    REFERENCES `mydb`.`Features` (`ID`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION);


-- -----------------------------------------------------
-- Table `mydb`.`Pubmed_Info`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Pubmed_Info` (
  `ID_pubmed` VARCHAR(50) NOT NULL,
  `Title` VARCHAR(200) NOT NULL,
  `Abstract` VARCHAR(2000) NOT NULL,
  `Pubmed_Infocol` VARCHAR(50) NOT NULL,
  `Authors` VARCHAR(100) NOT NULL,
  PRIMARY KEY (`ID_pubmed`, `Pubmed_Infocol`, `Authors`),
  CONSTRAINT `Pubmed-references`
    FOREIGN KEY (`ID_pubmed`)
    REFERENCES `mydb`.`References` (`Pubmed`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION);


-- -----------------------------------------------------
-- Table `mydb`.`Authors`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `mydb`.`Authors` (
  `ID_authors` VARCHAR(100) NOT NULL,
  `Name` VARCHAR(200) NOT NULL,
  `Affiliation` VARCHAR(200) NOT NULL,
  PRIMARY KEY (`ID_authors`),
  CONSTRAINT `Autores- pubmed`
    FOREIGN KEY (`ID_authors`)
    REFERENCES `mydb`.`Pubmed_Info` (`Authors`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION);


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
