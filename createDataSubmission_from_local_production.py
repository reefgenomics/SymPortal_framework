import os
import itertools
import subprocess
from dbApp.models import symportal_framework, data_set, reference_sequence, data_set_sample_sequence, analysis_type, analysis_group, data_set_sample, data_analysis, clade_collection, clade_collection_type
from multiprocessing import Queue, Process, Manager
from django import db
import pickle
import csv
import numpy as np
from collections import defaultdict
import shutil
import re
import json
import glob
from datetime import datetime

###### Generic functions ######
def readDefinedFileToList(filename):
    tempList = []
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite)-1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite)-1:
                writer.write(listToWrite[i])
            i += 1

def readByteObjectFromDefinedDirectory(directory):
    f = open(directory, 'rb')
    return pickle.load(f)

def writeByteObjectToDefinedDirectory(directory,object):
    f = open(directory , 'wb+')
    pickle.dump(object, f)

def createNoSpaceFastaFile(fastaList):
    tempList = []
    i = 0
    while i < len(fastaList):
        tempList.extend([fastaList[i].split('\t')[0], fastaList[i+1]])
        i += 2
    return tempList

def createDictFromFasta(fastaList):
    tempDict = {}
    i = 0
    while i < len(fastaList):
        sequence = fastaList[i][1:]
        tempDict[sequence] = fastaList[i+1]
        i += 2
    return tempDict

def createNewFile(pathtofile):
    try:
        os.makedirs(os.path.dirname(pathtofile))
    except FileExistsError:
        pass

    open(pathtofile, mode='w')

def writeLinesToFile(pathToFile, listoflines):
    with open(pathToFile, mode='a') as f:
        for line in listoflines:
            f.write(line + '\n')

def checkIfFileEmpty(filepath):
    count = 0
    with open(filepath, mode='r') as reader:
        for line in reader:
            if count != 0:
                return False  # not empty
            count += 1
    return True
###############################

def logQCErrorAndContinue(datasetsampleinstanceinq, samplename, errorreason):
    print('Error in processing sample: {}'.format(samplename))
    datasetsampleinstanceinq.finalUniqueSeqNum = 0
    datasetsampleinstanceinq.finalTotSeqNum = 0
    datasetsampleinstanceinq.initialProcessingComplete = True
    datasetsampleinstanceinq.errorInProcessing = True
    datasetsampleinstanceinq.errorReason = errorreason
    datasetsampleinstanceinq.save()

    return

def worker(input, output, wkd, dataSubID, e_val_collection_dict, reference_db_name):
    dataSubInQ = data_set.objects.get(id=dataSubID)
    for contigPair in iter(input.get, 'STOP'):
        sampleName = contigPair.split('\t')[0].replace('[dS]','-')

        dataSetSampleInstanceInQ = data_set_sample.objects.get(name=sampleName, dataSubmissionFrom=dataSubInQ)
        # Only process samples that have not already had this done.
        # This should be handy in the event of crashed midprocessing
        if dataSetSampleInstanceInQ.initialProcessingComplete == False:
            ###### NB We will always crop with the SYMVAR primers as they produce the shortest product
            primerFwdSeq = 'GAATTGCAGAACTCCGTGAACC'  # Written 5'-->3'
            primerRevSeq = 'CGGGTTCWCTTGTYTGACTTCATGC'  # Written 5'-->3'

            oligoFile = [
                r'#SYM_VAR_5.8S2',
                'forward\t{0}'.format(primerFwdSeq),
                r'#SYM_VAR_REV',
                'reverse\t{0}'.format(primerRevSeq)
            ]


            #Initial Mothur QC, making contigs, screening for ambiguous calls and homopolymers
            # Uniqueing, discarding <2 abundance seqs, removing primers and adapters

            print('Sample: {0}'.format(sampleName))
            currentDir = r'{0}/{1}/'.format(wkd, sampleName)
            os.makedirs(currentDir, exist_ok=True)
            stabilityFile = [contigPair]
            stabilityFileName = r'{0}{1}'.format(sampleName,'stability.files')
            rootName = r'{0}stability'.format(sampleName)
            stabilityFilePath = r'{0}{1}'.format(currentDir,stabilityFileName)
            writeListToDestination(stabilityFilePath, stabilityFile)
            # Write oligos file to directory
            writeListToDestination('{0}{1}'.format(currentDir, 'primers.oligos'), oligoFile)
            # NB mothur is working very strangely with the python subprocess command. For some
            # reason it is adding in an extra 'mothur' before the filename in the input directory
            # reason it is adding in an extra 'mothur' before the filename in the input directory
            # As such we will have to enter all of the paths to files absolutely

            # I am going to have to implement a check here that looks to see if the sequences are reverse complement or not.
            # Mothur pcr.seqs does not check to see if this is a problem. You simply get all of your seqs thrown out


            mBatchFile = [
                r'set.dir(input={0})'.format(currentDir),
                r'set.dir(output={0})'.format(currentDir),
                r'make.contigs(file={}{})'.format(currentDir, stabilityFileName),
                r'summary.seqs(fasta={}{}.trim.contigs.fasta)'.format(currentDir, rootName),
                r'screen.seqs(fasta={0}{1}.trim.contigs.fasta, group={0}{1}.contigs.groups, maxambig=0, maxhomop=5)'.format(
                    currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.fasta)'.format(currentDir, rootName),
                r'unique.seqs(fasta={0}{1}.trim.contigs.good.fasta)'.format(currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.fasta, name={0}{1}.trim.contigs.good.names)'.format(
                    currentDir, rootName),
                r'split.abund(cutoff=2, fasta={0}{1}.trim.contigs.good.unique.fasta, name={0}{1}.trim.contigs.good.names, group={0}{1}.contigs.good.groups)'.format(
                    currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta, name={0}{1}.trim.contigs.good.abund.names)'.format(
                    currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.rare.fasta, name={0}{1}.trim.contigs.good.rare.names)'.format(
                    currentDir, rootName),
                r'pcr.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta, group={0}{1}.contigs.good.abund.groups, name={0}{1}.trim.contigs.good.abund.names, oligos={0}primers.oligos, pdiffs=2, rdiffs=2)'.format(
                    currentDir, rootName)
            ]

            mBatchFilePath = r'{0}{1}{2}'.format(currentDir, 'mBatchFile', sampleName)
            writeListToDestination(mBatchFilePath, mBatchFile)

            # TODO Alejandros samples are causing us problems as they are presenting lots of errors when they
            # are running. We will try to intercept the error messages and kill the process if we get them
            #TODO if we get the blank fasta name error that causes Mothur to stall then we can take the contig
            # pair out and run it through pandaseq to get the contigs and then put it back into mothur.
            # This could be the way we save Alejandros 120 samples that seem to be causing an issue.
            error = False

            with subprocess.Popen(['mothur', '{0}'.format(mBatchFilePath)], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    # print(line)
                    if '[WARNING]: Blank fasta name, ignoring read.' in line:

                        #DEBUGing

                        p.terminate()
                        errorReason = 'Blank fasta name'
                        logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                        error = True
                        output.put(sampleName)
                        break
            apples = 'pears'
                # # It may be that we make it through the above mothur run without receiving the specific error message
                # # but an error may still have occoured so we should leave the following line in here.
                # if not error and p.returncode == 1:
                #     errorReason = 'make.contig but not Blank fasta name'
                #     logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                #     output.put(sampleName)
                #     continue
            if error:
                continue
            apples = 'pears'



            # Here check the outputted files to see if they are reverse complement or not by running the pcr.seqs and checking the results

            # Check to see if there are sequences in the PCR output file
            lastSummary = readDefinedFileToList(
                '{}{}.trim.contigs.good.unique.abund.pcr.fasta'.format(currentDir, rootName))
            if len(lastSummary) == 0:  # If this file is empty
                # Then these sequences may well be reverse complement so we need to try to rev first
                mBatchRev = [
                    r'reverse.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.fasta)'.format(currentDir, rootName),
                    r'pcr.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.rc.fasta, group={0}{1}.contigs.good.abund.groups, name={0}{1}.trim.contigs.good.abund.names, oligos={0}primers.oligos, pdiffs=2, rdiffs=2)'.format(
                        currentDir, rootName)
                ]
                mBatchFilePath = r'{0}{1}{2}'.format(currentDir, 'mBatchFile', sampleName)
                writeListToDestination(mBatchFilePath, mBatchRev)
                completedProcess = subprocess.run(
                    ['mothur', r'{0}'.format(mBatchFilePath)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                # At this point the sequences will be reversed and they will have been renamed so we
                # can just change the name of the .rc file to the orignal .fasta file that we inputted with
                # This way we don't need to change the rest of the mothur pipe.
                subprocess.run([r'mv', r'{0}{1}.trim.contigs.good.unique.abund.rc.pcr.fasta'.format(currentDir,rootName), r'{0}{1}.trim.contigs.good.unique.abund.pcr.fasta'.format(currentDir,rootName)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            lastSummary = readDefinedFileToList(
                '{}{}.trim.contigs.good.unique.abund.pcr.fasta'.format(currentDir, rootName))
            if len(lastSummary) == 0:  # If this file is still empty, then the problem was not solved by reverse complementing


                errorReason = 'error in inital QC'
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                output.put(sampleName)
                continue


            mBatchFileContinued = [
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.fasta, name={0}{1}.trim.contigs.good.abund.pcr.names)'.format(
                    currentDir, rootName),
                r'unique.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.fasta, name={0}{1}.trim.contigs.good.abund.pcr.names)'.format(
                    currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.unique.fasta, name={0}{1}.trim.contigs.good.unique.abund.pcr.names)'.format(
                    currentDir, rootName)
            ]

            mBatchFilePath = r'{0}{1}{2}'.format(currentDir, 'mBatchFile', sampleName)
            writeListToDestination(mBatchFilePath, mBatchFileContinued)
            completedProcess = subprocess.run(
                ['mothur', r'{0}'.format(mBatchFilePath)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #TODO examine which name file it is that we will need to break down in order to get the value for the
            # absolute number of sequences after the QC.

            if completedProcess.returncode == 1:
                errorReason = 'error in inital QC'
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                output.put(sampleName)
                continue
            print('{}: initial mothur complete'.format(sampleName))

            # Check to see if there are sequences in the PCR output file
            try:
                lastSummary = readDefinedFileToList(
                    '{}{}.trim.contigs.good.unique.abund.pcr.unique.fasta'.format(currentDir, rootName))
                if len(lastSummary) == 0:  # If this file is empty
                    errorReason = 'error in inital QC'
                    logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                    output.put(sampleName)
                    continue
            except: # If there is no file then we can assume sample has a problem
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName)
                continue


            # Get number of sequences after make.contig
            lastSummary = readDefinedFileToList('{}{}.trim.contigs.summary'.format(currentDir, rootName))
            number_of_seqs_contig_absolute = len(lastSummary) - 1
            dataSetSampleInstanceInQ.initialTotSeqNum = number_of_seqs_contig_absolute

            # Get number of sequences after unique
            lastSummary = readDefinedFileToList('{}{}.trim.contigs.good.unique.abund.pcr.unique.summary'.format(currentDir, rootName))
            number_of_seqs_contig_unique = len(lastSummary) - 1
            dataSetSampleInstanceInQ.initialUniqueSeqNum = number_of_seqs_contig_unique

            # Get absolute number of sequences after after sequence QC
            last_summary = readDefinedFileToList('{}{}.trim.contigs.good.unique.abund.pcr.unique.summary'.format(currentDir, rootName))
            absolute_count = 0
            for line in last_summary[1:]:
                absolute_count += int(line.split('\t')[6])
            dataSetSampleInstanceInQ.post_seq_qc_absolute_num_seqs = absolute_count
            dataSetSampleInstanceInQ.save()

            if sampleName == 'P7-F05_P7-F05_N705-S520':
                apples = 'asdf'

            print('Initial mothur complete')
            # Each sampleDataDir should contain a set of .fasta, .name and .group files that we can use to do local blasts with

            ncbircFile = []
            db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))
            ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])





            # Run local blast of all seqs and determine clade. Discard seqs below evalue cutoff and write out new fasta, name, group and clade dict
            print('Verifying seqs are Symbiodinium and determining clade: {}'.format(rootName.replace('stability','')))

            #write the .ncbirc file that gives the location of the db
            writeListToDestination("{0}.ncbirc".format(currentDir), ncbircFile)

            #Read in the fasta, name and group files and convert to dics
            fastaFile = readDefinedFileToList('{0}{1}.trim.contigs.good.unique.abund.pcr.unique.fasta'.format(currentDir, rootName))
            uniqueFastaFile = createNoSpaceFastaFile(fastaFile)
            writeListToDestination('{}blastInputFasta.fa'.format(currentDir), uniqueFastaFile)
            fastaDict = createDictFromFasta(uniqueFastaFile)
            nameFile = readDefinedFileToList('{0}{1}.trim.contigs.good.unique.abund.pcr.names'.format(currentDir, rootName))
            nameDict = {a.split('\t')[0]: a for a in nameFile}

            groupFile = readDefinedFileToList('{0}{1}.contigs.good.abund.pcr.groups'.format(currentDir, rootName))

            # Set up environment for running local blast

            blastOutputPath = r'{}blast.out'.format(currentDir)
            outputFmt = "6 qseqid sseqid staxids evalue"
            inputPath = r'{}blastInputFasta.fa'.format(currentDir)
            os.chdir(currentDir)

            # Run local blast
            # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
            completedProcess = subprocess.run(['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', reference_db_name,
                 '-max_target_seqs', '1', '-num_threads', '1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('{}: blast complete'.format(sampleName))
            # Read in blast output
            blastOutputFile = readDefinedFileToList(r'{}blast.out'.format(currentDir))
            blastDict = {a.split('\t')[0]: a.split('\t')[1][-1] for a in blastOutputFile}
            throwAwaySeqs = []
            ###### Uncomment for blasting QC
            #######blastedSeqs = []
            #####NB it turns out that the blast will not always return a match.
            # If a match is not returned it means the sequence did not have a significant match to the seqs in the db
            #Add any seqs that did not return a blast match to the throwAwaySeq list
            diff = set(fastaDict.keys()) - set(blastDict.keys())
            throwAwaySeqs.extend(list(diff))

            ## 030518 We are starting to throw away Symbiodinium sequences here, especially in the non-coral samples
            # I think we will need to severely relax the e value cut off in order to incorporate more sequences

            # NB note that the blast results sometime return several matches for the same seq.
            # as such we will use the already_processed_blast_seq_resulst to make sure that we only
            # process each sequence once.
            already_processed_blast_seq_result = []
            for line in blastOutputFile:
                seqInQ = line.split('\t')[0]
                if seqInQ in already_processed_blast_seq_result:
                    continue
                already_processed_blast_seq_result.append(seqInQ)
                try:
                    evaluePower = int(line.split('\t')[3].split('-')[1])
                    if evaluePower < 100:  # evalue cut off, collect sequences that don't make the cut
                        throwAwaySeqs.append(seqInQ)
                        # incorporate the size cutoff here that would normally happen below
                        if 184 < len(fastaDict[seqInQ]) < 310:
                            if fastaDict[seqInQ] in e_val_collection_dict.keys():
                                e_val_collection_dict[fastaDict[seqInQ]] += 1
                            else:
                                e_val_collection_dict[fastaDict[seqInQ]] = 1
                except:
                    throwAwaySeqs.append(seqInQ)
                    if 184 < len(fastaDict[seqInQ]) < 310:
                        if fastaDict[seqInQ] in e_val_collection_dict.keys():
                            e_val_collection_dict[fastaDict[seqInQ]] += 1
                        else:
                            e_val_collection_dict[fastaDict[seqInQ]] = 1

            #have a look at the above code and see how we can get a total count for all of the sequences
            # that were thrown away. I think we'll likely have to look for a name file to be able to look
            # the thrown away sequences up in.
            # Add number of absolute sesquences that are not Symbiodinium

            # NB it turns out that sometimes a sequence is returned in the blast results twice! This was messing up
            # our meta-analysis reporting. This will be fixed by working with sets of the throwaway sequences
            temp_count = 0
            for seq_name in list(set(throwAwaySeqs)):
                temp_count += len(nameDict[seq_name].split('\t')[1].split(','))
            dataSetSampleInstanceInQ.non_sym_absolute_num_seqs = temp_count
            # Add details of non-symbiodinium unique seqs
            dataSetSampleInstanceInQ.nonSymSeqsNum = len(set(throwAwaySeqs))
            dataSetSampleInstanceInQ.save()

            ###### Uncomment for blasting QC
            ########notBlasted = list(set(blastedSeqs)-set(nameDict.keys()))
            ########print('notblasted: ' + str(len(notBlasted)))

            # Output new fasta, name and group files that don't contain seqs that didn't make the cut

            print('Discarding {} unique sequences for evalue cutoff violations'.format(str(len(throwAwaySeqs))))
            newFasta = []
            newName = []
            newGroup = []
            cladalDict = {}
            count = 0
            listOfBadSeqs = []
            for line in groupFile:
                sequence = line.split('\t')[0]
                if sequence not in throwAwaySeqs:
                    newGroup.append(line)
                    # The fastaDict is only meant to have the unique seqs in so this will go to 'except' a lot. This is OK and normal
                    try:
                        newFasta.extend(['>{}'.format(sequence), fastaDict[sequence]])
                    except:
                        pass

                    try:
                        newName.append(nameDict[sequence])
                        listOfSameSeqNames = nameDict[sequence].split('\t')[1].split(',')
                        clade = blastDict[sequence]

                        for seqName in listOfSameSeqNames:
                            cladalDict[seqName] = clade
                    except:
                        pass
            # Now write the files out
            if not newFasta:
                # Then the fasta is blank and we have got no good Symbiodinium seqs
                errorReason = 'No Symbiodinium sequences left after blast annotation'
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                output.put(sampleName)
                continue
            print('{}: non-Symbiodinium sequences binned'.format(sampleName))
            writeListToDestination('{0}{1}.trim.contigs.good.unique.abund.pcr.blast.fasta'.format(currentDir, rootName), newFasta)
            writeListToDestination('{0}{1}.trim.contigs.good.abund.pcr.blast.names'.format(currentDir, rootName), newName)
            writeListToDestination('{0}{1}.contigs.good.abund.pcr.blast.groups'.format(currentDir, rootName), newGroup)
            writeByteObjectToDefinedDirectory('{0}{1}.cladeDict.dict'.format(currentDir, rootName), cladalDict)
            # At this point we have the newFasta, newName, newGroup. These all only contain sequences that were
            # above the blast evalue cutoff.
            # We also have the cladalDict for all of these sequences.summary

            # Now finish off the mothur analyses by discarding by size range
            # Have to find how big the average seq was and then discard 50 bigger or smaller than this
            # Read in last summary file
            #TODO check to see if this is the correct last Summary file
            # TODO maybe consider doing absolute size range cutoff
            # I am now going to switch this to an absolute size range as I am having problems with Mani's sequences.
            # For some reason he is having an extraordinarily high number of very short sequence (i.e. 15bp long).
            # These are not being thrown out in the blast work. As such the average is being thrown off. and means that our
            # upper size limit is only about 200.
            # I have calculated the averages of each of the clades for our reference sequences so far
            '''Clade A 234.09815950920245
                Clade B 266.79896907216494
                Clade C 261.86832986832985
                Clade D 260.44158075601376
                 '''
            # I will take our absolute cutoffs from these numbers (+- 50 bp) so 184-310
            secondmBatchFilePathList = []
            lastSummary = readDefinedFileToList('{0}{1}.trim.contigs.good.unique.abund.pcr.unique.summary'.format(currentDir, rootName))
            sum = 0
            for line in lastSummary[1:]:
                sum += int(line.split('\t')[3])
            average = int(sum/len(lastSummary))
            cutOffLower = 184
            cutOffUpper = 310


            secondmBatchFile = [
                r'set.dir(input={0})'.format(currentDir),
                r'set.dir(output={0})'.format(currentDir),
                r'screen.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.blast.fasta, name={0}{1}.trim.contigs.good.abund.pcr.blast.names, group={0}{1}.contigs.good.abund.pcr.blast.groups,  minlength={2}, maxlength={3})'.format(currentDir,rootName,cutOffLower,cutOffUpper),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.fasta, name={0}{1}.trim.contigs.good.abund.pcr.blast.good.names)'.format(currentDir,rootName),
                r'unique.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.fasta, name={0}{1}.trim.contigs.good.abund.pcr.blast.good.names)'.format(
                    currentDir, rootName),
                r'summary.seqs(fasta={0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.unique.fasta, name={0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.names)'.format(
                    currentDir, rootName),
            ]
            mBatchFilePath = r'{0}{1}{2}'.format(currentDir, 'mBatchFile', rootName)
            secondmBatchFilePathList.append(mBatchFilePath)

            writeListToDestination(mBatchFilePath, secondmBatchFile)
            completedProcess = subprocess.run(['mothur', r'{0}'.format(mBatchFilePath)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if completedProcess.returncode == 1:
                errorReason = 'No Symbiodinium sequences left after size screening'
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName, errorReason)
                output.put(sampleName)
                continue


            #### Here make cladally separated fastas

            try:
                fastaFile = readDefinedFileToList('{0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.unique.fasta'.format(currentDir, rootName))
                nameFile = readDefinedFileToList('{0}{1}.trim.contigs.good.unique.abund.pcr.blast.good.names'.format(currentDir, rootName))
            except:
                logQCErrorAndContinue(dataSetSampleInstanceInQ, sampleName)
                continue

            print('{}: final Mothur completed'.format(sampleName))
            fastaDict = createDictFromFasta(fastaFile)

            nameDict = {a.split('\t')[0]: a for a in nameFile}
            cladeDict = readByteObjectFromDefinedDirectory('{0}{1}.cladeDict.dict'.format(currentDir, rootName))
            cladeDirs = []
            cladeFastas = {}
            for line in nameFile:
                sequence = line.split('\t')[0]
                clade = cladeDict[sequence]
                if clade in cladeDirs:#Already have a dir for it
                    cladeFastas[clade][0].extend(['>{}'.format(sequence), fastaDict[sequence]])
                    cladeFastas[clade][1].append(nameDict[sequence])

                else: #Make dir and add empty fasta list to cladeFastas
                    cladeFastas[clade] = ([],[])
                    cladeDirs.append(clade)
                    os.makedirs(r'{}{}'.format(currentDir,clade), exist_ok=True)
                    cladeFastas[clade][0].extend(['>{}'.format(sequence), fastaDict[sequence]])
                    cladeFastas[clade][1].append(nameDict[sequence])

            total_debug_absolute = 0
            total_debug_unique = 0
            for someclade in cladeDirs:
                ### Debug###

                # work out the absolute number of sequences and unique sequences from these file and compare to the below
                temp_name = cladeFastas[someclade][1]
                total_debug_unique += len(temp_name)
                for temp_line in temp_name:
                    total_debug_absolute += len(temp_line.split('\t')[1].split(','))
                ####
                writeListToDestination(r'{0}{1}/{2}.QCed.clade{1}.fasta'.format(currentDir,someclade,rootName.replace('stability','')),cladeFastas[someclade][0])
                writeListToDestination(r'{0}{1}/{2}.QCed.clade{1}.names'.format(currentDir, someclade, rootName.replace('stability','')), cladeFastas[someclade][1])

            # Here we have cladaly sorted fasta and name file in new directory

            # now populate the data set sample with the qc meta-data
            # get unique seqs remaining
            dataSetSampleInstanceInQ.finalUniqueSeqNum = len(nameDict)
            #Get total number of sequences
            count = 0
            for nameKey in nameDict.keys():
                count += len(nameDict[nameKey].split('\t')[1].split(','))
            dataSetSampleInstanceInQ.finalTotSeqNum = count
            # now get the seqs lost through size violations through subtraction
            dataSetSampleInstanceInQ.size_violation_absolute = dataSetSampleInstanceInQ.post_seq_qc_absolute_num_seqs - dataSetSampleInstanceInQ.finalTotSeqNum - dataSetSampleInstanceInQ.non_sym_absolute_num_seqs
            dataSetSampleInstanceInQ.size_violation_unique = dataSetSampleInstanceInQ.initialUniqueSeqNum - dataSetSampleInstanceInQ.finalUniqueSeqNum - dataSetSampleInstanceInQ.nonSymSeqsNum

            # Now update the data_set_sample instance to set initialProcessingComplete to True
            dataSetSampleInstanceInQ.initialProcessingComplete = True
            dataSetSampleInstanceInQ.save()
            print('{}: initial processing complete'.format(sampleName))

            os.chdir(currentDir)
            fileList = [f for f in os.listdir(currentDir) if f.endswith((".names", ".fasta", ".qual", ".summary", ".oligos",
                                                             ".accnos", ".files", ".groups", ".logfile", ".dict", ".fa",
                                                             ".out"))]
            for f in fileList:
                os.remove(f)

    return

def deuniqueFiles(wkd, ID, numProc):

    # Create mothur batch for each .fasta .name pair to be deuniqued
    # Put in directory list, run via multiprocessing
    samplesCollection = data_set_sample.objects.filter(dataSubmissionFrom=data_set.objects.get(id=ID))
    mBatchFilePathList = []
    for dataSetSampleInstance in samplesCollection: # For each samples directory
        sampleName = dataSetSampleInstance.name
        fullPath = '{}/{}'.format(wkd, sampleName)

        #http: // stackoverflow.com / questions / 3207219 / how - to - list - all - files - of - a - directory
        listOfDirs = []
        for (dirpath, dirnames, filenames) in os.walk(fullPath):
            listOfDirs.extend(dirnames)
            break
        for directory in listOfDirs:# for each cladal directory
            fastaFilePath = ''
            nameFilePath = ''
            pathToDir = '{0}/{1}'.format(fullPath, directory)
            cladeName = directory
            # For each of the files in each of the Cladal directories
            listOfFiles = []
            for (dirpath, dirnames, filenames) in os.walk(pathToDir):
                listOfFiles.extend(filenames)
                break

            for files in listOfFiles:
                if '.fasta' in files and '.redundant' not in files:
                    fastaFilePath = '{0}/{1}'.format(pathToDir, files)
                elif '.names' in files:
                    nameFilePath = '{0}/{1}'.format(pathToDir, files)

            # Build a quick mBatchFile
            mBatchFile = [
                r'set.dir(input={0}/)'.format(pathToDir),
                r'set.dir(output={0}/)'.format(pathToDir),
                r'deunique.seqs(fasta={0}, name={1})'.format(fastaFilePath, nameFilePath)
            ]
            mBatchFilePath = '{0}/{1}'.format(pathToDir, '{0}.{1}.{2}'.format(sampleName, cladeName, 'mBatchFile'))
            writeListToDestination(mBatchFilePath, mBatchFile)
            mBatchFilePathList.append(mBatchFilePath)

    # Create the queues that will hold the mBatchFile paths
    taskQueue = Queue()
    doneQueue = Queue()

    for mBatchFilePath in mBatchFilePathList:
        taskQueue.put(mBatchFilePath)


    for n in range(numProc):
        taskQueue.put('STOP')

    allProcesses = []

    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    for n in range(numProc):
        p = Process(target=deuniqueWorker, args=(taskQueue, doneQueue))
        allProcesses.append(p)
        p.start()

    # Collect the list of deuniqued directories to use for MED analyses
    listOfDeuniquedFastaPaths = []
    for i in range(len(mBatchFilePathList)):
        listOfDeuniquedFastaPaths.append(doneQueue.get())

    for p in allProcesses:
        p.join()

    ### DEBUG ###

    #############

    return listOfDeuniquedFastaPaths

def deuniqueWorker(input, output):

    # This currently works through a list of paths to batch files to be uniques.
    # But at each of these locations once the modified deuniqued file has been written we can then perform the MED
    # analysis on the file in each of the directories.
    # We also want to be able to read in the results of the MED but we will not be able to do that as MP so we
    # will have to save the list of directories and go through them one by one to create the sequences

    for mBatchFilePath in iter(input.get, 'STOP'):
        found = True

        # Run the dunique
        completedProcess = subprocess.run(['mothur', r'{0}'.format(mBatchFilePath)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Modify the deuniqued fasta to append sample name to all of the sequences
        # Get list of files in directory
        deuniquedFasta = []
        cwd = os.path.dirname(mBatchFilePath)
        sampleName = cwd.split('/')[-2]
        # Replace '_' in name as MED uses text up to first underscore as sample name
        # This shouldn't be necessary
        #sampleName = sampleName.replace('_', '-')
        listOfFiles = []
        for (dirpath, dirnames, filenames) in os.walk(cwd):
            listOfFiles.extend(filenames)
            break
        pathToFile = None
        for file in listOfFiles:
            if '.redundant' in file: # Then this is the deuniqued fasta
                pathToFile = '{0}/{1}'.format(cwd, file)

                break
        deuniquedFasta = readDefinedFileToList(pathToFile)
        deuniquedFasta = ['{0}{1}_{2}'.format(a[0],sampleName,a[1:].replace('_','-')) if a[0] == '>' else a for a in deuniquedFasta]
        #write the modified deuniquedFasta to list
        writeListToDestination(pathToFile, deuniquedFasta)
        # Put the path to the deuniqued fasta into the output list for use in MED analyses
        output.put('{}/{}/'.format(os.path.dirname(pathToFile), 'MEDOUT'))

        try:
            # The fasta that we want to pad and MED is the 'file'
            completedProcess = subprocess.run([r'o-pad-with-gaps', r'{}'.format(pathToFile)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if completedProcess.returncode == 1:
                pear = 'appe;'
            # Now run MED
            listOfFiles = []
            for (dirpath, dirnames, filenames) in os.walk(cwd):
                listOfFiles.extend(filenames)
                break
            for file in listOfFiles:
                if 'PADDED' in file:
                    pathToFile = '{0}/{1}'.format(cwd, file)
                    break
            MEDOutDir = '{}/{}/'.format(cwd, 'MEDOUT')
            os.makedirs(MEDOutDir, exist_ok=True)
            completedProcess = subprocess.run(
                [r'decompose', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o',
                 MEDOutDir, pathToFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if completedProcess.returncode == 1:
                pear = 'appe;'
        except:
            breakPpoint = 'point'

def checkIfSeqInQHadRefSeqMatch(seqInQ, nodeName, refSeqIdDict, nodeToRefDict):
    # We use this to look to see if there is an equivalent refSeq Sequence for the sequence in question
    # This take into account whether the seqInQ could be a subset or super set of one of the
    # refSeq.sequences
    # Will return false if no refSeq match is found

    # first check to see if seq is found
    if seqInQ in refSeqIdDict:  # Found actual seq in dict
        nodeToRefDict[nodeName] = refSeqIdDict[seqInQ]
        print('Assigning existing reference sequence {} to MED node {}'.format(reference_sequence.objects.get(id=refSeqIdDict[seqInQ]).name, nodeName))
        return True
    elif 'A' + seqInQ in refSeqIdDict:  # This was a seq shorter than refseq but we can associate it to this ref seq
        nodeToRefDict[nodeName] = refSeqIdDict['A' + seqInQ]
        print('Assigning existing reference sequence {} to MED node {}'.format(reference_sequence.objects.get(id=refSeqIdDict['A' + seqInQ]).name, nodeName))
        return True
    else:  # This checks if either the seq in question is found in the sequence of a reference_sequence
        # or if the seq in question is bigger than a refseq sequence and is a super set of it
        # In either of these cases we should consider this a match and use the refseq matched to.
        # This might be very coputationally expensive but lets give it a go

        for seqKey in refSeqIdDict.keys():
            if seqInQ in seqKey or seqKey in seqInQ:
                # Then this is a match
                nodeToRefDict[nodeName] = refSeqIdDict[seqKey]
                print('Assigning existing reference sequence {} to MED node {}'.format(reference_sequence.objects.get(id=refSeqIdDict[seqKey]).name, nodeName))
                return True
    return False

def collateFilesForMed(listofdeuniquedfastapaths, wkd):
    # To avoid a memory crash we append each deuniqued fasta directly to the master fasta on disk
    # This way we don't hold the master fasta in memory as a list
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    # wkd = '/'.join(listofdeuniquedfastapaths[0].split('/')[:5])
    for clade in cladeList:
        createNewFile('{}/deuniquedFastaForMED/redundant.{}.fasta'.format(wkd, clade))
    # list will hold all of the fastas
    #cladeSepFastaList = [[] for a in cladeList]
    for deUFastaPath in listofdeuniquedfastapaths:
        deuniquedFasta = readDefinedFileToList(deUFastaPath)
        clade = deUFastaPath.split('/')[-2]
        writeLinesToFile('{}/deuniquedFastaForMED/redundant.{}.fasta'.format(wkd, clade), deuniquedFasta)

    # Count number of seqs in each file
    # If empty, then delete
    for clade in cladeList:
        fastaPath = '{}/deuniquedFastaForMED/redundant.{}.fasta'.format(wkd, clade)
        if checkIfFileEmpty(fastaPath): # Delete the MED input clade files that are empty
            os.remove(fastaPath)

    return

def runMED(wkd):
    # For each of the files in the wkd/deuniquedFastaForMed directory
    # Pad gaps
    # Decompose


    # Get files
    # and pad

    listOfFiles = []
    for (dirpath, dirnames, filenames) in os.walk('{}/deuniquedFastaForMED/'.format(wkd)):
        listOfFiles.extend(filenames)
        break
    for file in listOfFiles:
        if '.fasta' in file and 'PADDED' not in file:
            pathToFile = '{}/deuniquedFastaForMED/{}'.format(wkd, file)
            # Sanity check at debug to how sample information is extracted
            #completedProcess = subprocess.run([r'o-get-sample-info-from-fasta', r'{0}'.format(pathToFile)])
            completedProcess = subprocess.run([r'o-pad-with-gaps', r'{0}'.format(pathToFile)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #Decompose
    listOfFiles = []
    for (dirpath, dirnames, filenames) in os.walk('{}/deuniquedFastaForMED/'.format(wkd)):
        listOfFiles.extend(filenames)
        break
    for file in listOfFiles:
        if 'PADDED' in file: # Then this is one of the padded fastas
            clade = file.split('.')[1]
            pathToFile = '{}/deuniquedFastaForMED/{}'.format(wkd, file)
            outputDir = '{}/deuniquedFastaForMED/{}'.format(wkd, clade)
            os.makedirs(outputDir, exist_ok=True)
            completedProcess = subprocess.run([r'decompose', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o', outputDir, pathToFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return

def processMEDDataDirectCCDefinition(wkd, ID, MEDDirs):
    # We are going to change this so that we go to each of the MEDDirs, which represent the clades within samples
    # that have had MED analyses run in them and we are going to use the below code to populate sequences to
    # the CCs and samples
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dir in MEDDirs: # For each of the directories where we did MED
        os.chdir(dir)
        # Get the sample
        sampleName = dir.split('/')[-4]
        sampleObject = data_set_sample.objects.get(dataSubmissionFrom=data_set.objects.get(id=ID), name=sampleName)
        # Get the clade
        clade = dir.split('/')[-3]
        print('Populating {} with clade {} sequences'.format(sampleName, clade))
        # Read in the node file
        try:
            nodeFile = readDefinedFileToList('NODE-REPRESENTATIVES.fasta')
        except:
            continue
        nodeFile = [line.replace('-', '') if line[0] != '>' else line for line in nodeFile]

        # Create nodeToRefDict that will be populated
        nodeToRefDict = {}

        # This is a dict of key = reference_sequence.sequence value = reference_sequence.id for all refseqs
        # We will use this to see if the sequence in question has a match , or is found in (this is key
        # as some of the seqs are one bp smaller than the reference seqs) there reference sequences
        refSeqDict = {refSeq.sequence: refSeq.id for refSeq in reference_sequence.objects.all()}

        ############ ASSOCIATE MED NODES TO EXISITING REFSEQS OR CREATE NEW REFSEQS #########
        # Look up seq of each with reference_sequence table
        # See if the seq in Q matches a reference_sequence, if so, associate
        listOfRefSeqs = []
        for i in range(len(nodeFile)):
            # We were having a problem here where some of the seqs were 1bp shorter than the reference seqs
            # As such they werent matching to the refernceSequence object e.g. to C3 but when we do the
            # blast they come up with C3 as their clostest and perfect match
            # To fix this we will run checkIfSeqInQHadRefSeqMatch

            if nodeFile[i][0] == '>':  # Then this is a def line
                sequenceInQ = nodeFile[i + 1]
                nodeNameInQ = nodeFile[i][1:].split('|')[0]
                # If True then node name will already have been associated to nodeToRefDict
                # and no need to do anything else
                found = checkIfSeqInQHadRefSeqMatch(seqInQ=sequenceInQ, refSeqIdDict=refSeqDict, nodeName=nodeNameInQ,
                                                    nodeToRefDict=nodeToRefDict)

                if found == False:
                    ####listToBeBlasted.extend([nodeFile[i], sequenceInQ])
                    newreferenceSequence = reference_sequence(clade=clade, sequence=sequenceInQ)
                    newreferenceSequence.save()
                    listOfRefSeqs.append(newreferenceSequence)
                    refSeqDict[newreferenceSequence.sequence] = newreferenceSequence.id
                    nodeToRefDict[nodeNameInQ] = newreferenceSequence.id

                    print('Assigning new reference sequence {} to MED node {}'.format(newreferenceSequence.name,
                                                                                      nodeFile[i][1:].split('|')[0]))
        ########################################################################################

        # # Here we have a refSeq associated to each of the seqs found and we can now create dataSetSampleSequences that have associated referenceSequences
        # So at this point we have a reference_sequence associated with each of the nodes
        # Now it is time to define clade collections
        # Open count table as list of lists
        countArray = []
        nodes = []
        samples = []
        # this creates countArray which is a 2D list
        with open('MATRIX-COUNT.txt') as f:
            reader = csv.reader(f, delimiter='\t')
            countArray = list(reader)
        # get Nodes from first list
        nodes = countArray[0][1:]
        # del nodes
        del countArray[0]
        # get samples from first item of each list
        # del samples to leave raw numerical
        for i in range(len(countArray)):
            samples.append(countArray[i][0])
            del countArray[i][0]
        # convert to np array
        countArray = np.array(countArray)
        countArray = countArray.astype(np.int)
        # for each node in each sample create dataSetSampleSeq with foreign key to referenceSeq and data_set_sample
        # give it a foreign key to the reference Seq by looking up the seq in the dictionary made earlier and using the value to search for the referenceSeq



        for i in range(len(samples)):  # For each sample # There should only be one sample

            sampleObject = data_set_sample.objects.get(dataSubmissionFrom=data_set.objects.get(id=ID),
                                                     name=samples[i])
            # Add the metadata to the data_set_sample
            sampleObject.post_med_absolute += sum(countArray[i])
            sampleObject.post_med_unique += len(countArray[i])
            sampleObject.save()

            cladalSeqAbundanceCounter = [int(a) for a in json.loads(sampleObject.cladalSeqTotals)]
            # TODO This is where we need to tackle the issue of making sure we keep track of sequences in samples that
            # were not above the 200 threshold to be made into cladeCollections
            # We will simply add a list to the sampleObject that will be a sequence total for each of the clades
            # in order of cladeList
            if sum(countArray[i]) > 200:  # Then this sample has enough seqs in this clade to create a cladeCollection
                # Here we modify the cladalSeqTotals string of the sample object to add the sequence totals
                # for the given clade
                cladeIndex = cladeList.index(clade)
                tempInt = cladalSeqAbundanceCounter[cladeIndex]
                tempInt += sum(countArray[i])
                cladalSeqAbundanceCounter[cladeIndex] = tempInt
                sampleObject.cladalSeqTotals = json.dumps([str(a) for a in cladalSeqAbundanceCounter])
                sampleObject.save()


                dssList = []
                print('Populating sample {} with clade {} sequences'.format(samples[i], clade))
                # if > 200 then create the CC
                # don't
                newCC = clade_collection(clade=clade, dataSetSampleFrom=sampleObject)
                newCC.save()


                refSeqAbundanceCounter = defaultdict(int)
                for j in range(len(
                        nodes)):  # Only process this sample for this clade if there are enough seqs to make a cladeCollection
                    abundance = countArray[i][j]
                    if abundance > 0:
                        # I want to address a problem we are having here. Now that we have thorough checks to
                        # associate very similar sequences with indels by the primers to the same reference seq
                        # it means that multiple sequences within the same sample can have the same referenceseqs
                        # Due to the fact that we will in effect use the sequence of the reference seq rather
                        # than the dsss seq, we should consolidate all dsss seqs with the same reference seq
                        # so... we will create a counter that will keep track of the cumulative abundance associated with each reference_sequence
                        # and then create a dsss for each refSeq
                        refSeqAbundanceCounter[
                            reference_sequence.objects.get(id=nodeToRefDict[nodes[j]])] += abundance
                # > 200 associate a CC to the data_set_sample, else, don't
                # irrespective, associate a data_set_sample
                for refSeq in refSeqAbundanceCounter.keys():
                    dss = data_set_sample_sequence(referenceSequenceOf=refSeq,
                                                   cladeCollectionTwoFoundIn=newCC,
                                                   abundance=refSeqAbundanceCounter[refSeq])
                    dssList.append(dss)
                # Save all of the newly created dss
                data_set_sample_sequence.objects.bulk_create(dssList)
                # Get the ids of each of the dss and add create a string of them and store it as cc.footPrint
                # This way we can quickly get the footprint of the CC.
                # Sadly we can't get eh IDs from the list so we will need to re-query
                # Instead we add the ID of each refseq in the refSeqAbundanceCounter.keys() list
                newCC.footPrint = ','.join([str(refSeq.id) for refSeq in refSeqAbundanceCounter.keys()])
                newCC.save()
            # get rid of this.
            elif 0 < sum(countArray[i]) < 200:
                # Then these are sequences that will be outside of a clade Collection
                # For each sample object in the database we will have a list that will hold the abundances
                # Of sequences that are outside of cladeCollections

                # NB Bear in mind that we make a reference_sequence object for every node found irrespecitve
                # of if it is used in a cladeCollection or not which is great!
                # Latter on in the output we will be able to output all sequences that have been found so far.
                # Only sequence that are DIVs will have names assocciated to them.

                # Here we modify the cladalSeqTotals string of the sample object to add the sequence totals
                # for the given clade
                cladeIndex = cladeList.index(clade)
                tempInt = cladalSeqAbundanceCounter[cladeIndex]
                tempInt += sum(countArray[i])
                cladalSeqAbundanceCounter[cladeIndex] = tempInt
                sampleObject.cladalSeqTotals = json.dumps([str(a) for a in cladalSeqAbundanceCounter])
                sampleObject.save()

    return

def processMEDDataDirectCCDefinition_new_dss_structure(wkd, ID, MEDDirs):
    ''' Here we have modified the original method processMEDDataDirectCCDefinition'''
    # We are going to change this so that we go to each of the MEDDirs, which represent the clades within samples
    # that have had MED analyses run in them and we are going to use the below code to populate sequences to
    # the CCs and samples
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    for dir in MEDDirs: # For each of the directories where we did MED
        os.chdir(dir)
        # Get the sample
        sampleName = dir.split('/')[-4]
        data_set_sample_object = data_set_sample.objects.get(dataSubmissionFrom=data_set.objects.get(id=ID), name=sampleName)
        # Get the clade
        clade = dir.split('/')[-3]
        print('Populating {} with clade {} sequences'.format(sampleName, clade))
        # Read in the node file
        try:
            nodeFile = readDefinedFileToList('NODE-REPRESENTATIVES.fasta')
        except:
            continue
        nodeFile = [line.replace('-', '') if line[0] != '>' else line for line in nodeFile]

        # Create nodeToRefDict that will be populated
        nodeToRefDict = {}

        # This is a dict of key = reference_sequence.sequence value = reference_sequence.id for all refseqs
        # We will use this to see if the sequence in question has a match , or is found in (this is key
        # as some of the seqs are one bp smaller than the reference seqs) there reference sequences
        refSeqDict = {refSeq.sequence: refSeq.id for refSeq in reference_sequence.objects.all()}

        ############ ASSOCIATE MED NODES TO EXISITING REFSEQS OR CREATE NEW REFSEQS #########
        # Look up seq of each with reference_sequence table
        # See if the seq in Q matches a reference_sequence, if so, associate
        listOfRefSeqs = []
        for i in range(len(nodeFile)):
            # We were having a problem here where some of the seqs were 1bp shorter than the reference seqs
            # As such they werent matching to the refernceSequence object e.g. to C3 but when we do the
            # blast they come up with C3 as their clostest and perfect match
            # To fix this we will run checkIfSeqInQHadRefSeqMatch

            if nodeFile[i][0] == '>':  # Then this is a def line
                sequenceInQ = nodeFile[i + 1]
                nodeNameInQ = nodeFile[i][1:].split('|')[0]
                # If True then node name will already have been associated to nodeToRefDict
                # and no need to do anything else
                found = checkIfSeqInQHadRefSeqMatch(seqInQ=sequenceInQ, refSeqIdDict=refSeqDict, nodeName=nodeNameInQ,
                                                    nodeToRefDict=nodeToRefDict)

                if found == False:
                    ####listToBeBlasted.extend([nodeFile[i], sequenceInQ])
                    newreferenceSequence = reference_sequence(clade=clade, sequence=sequenceInQ)
                    newreferenceSequence.save()
                    listOfRefSeqs.append(newreferenceSequence)
                    refSeqDict[newreferenceSequence.sequence] = newreferenceSequence.id
                    nodeToRefDict[nodeNameInQ] = newreferenceSequence.id

                    print('Assigning new reference sequence {} to MED node {}'.format(newreferenceSequence.name,
                                                                                      nodeFile[i][1:].split('|')[0]))
        ########################################################################################

        # # Here we have a refSeq associated to each of the seqs found and we can now create dataSetSampleSequences that have associated referenceSequences
        # So at this point we have a reference_sequence associated with each of the nodes
        # Now it is time to define clade collections
        # Open count table as list of lists
        countArray = []
        nodes = []
        samples = []
        # this creates countArray which is a 2D list
        with open('MATRIX-COUNT.txt') as f:
            reader = csv.reader(f, delimiter='\t')
            countArray = list(reader)
        # get Nodes from first list
        nodes = countArray[0][1:]
        # del nodes
        del countArray[0]
        # get samples from first item of each list
        # del samples to leave raw numerical
        for i in range(len(countArray)):
            samples.append(countArray[i][0])
            del countArray[i][0]
        # convert to np array
        countArray = np.array(countArray)
        countArray = countArray.astype(np.int)
        # for each node in each sample create dataSetSampleSeq with foreign key to referenceSeq and data_set_sample
        # give it a foreign key to the reference Seq by looking up the seq in the dictionary made earlier and using the value to search for the referenceSeq



        for i in range(len(samples)):  # For each sample # There should only be one sample

            data_set_sample_object = data_set_sample.objects.get(dataSubmissionFrom=data_set.objects.get(id=ID),
                                                     name=samples[i])
            # Add the metadata to the data_set_sample
            data_set_sample_object.post_med_absolute += sum(countArray[i])
            data_set_sample_object.post_med_unique += len(countArray[i])
            data_set_sample_object.save()
            cladalSeqAbundanceCounter = [int(a) for a in json.loads(data_set_sample_object.cladalSeqTotals)]

            # This is where we need to tackle the issue of making sure we keep track of sequences in samples that
            # were not above the 200 threshold to be made into cladeCollections
            # We will simply add a list to the sampleObject that will be a sequence total for each of the clades
            # in order of cladeList

            # Here we modify the cladalSeqTotals string of the sample object to add the sequence totals
            # for the given clade
            cladeIndex = cladeList.index(clade)
            tempInt = cladalSeqAbundanceCounter[cladeIndex]
            tempInt += sum(countArray[i])
            cladalSeqAbundanceCounter[cladeIndex] = tempInt
            data_set_sample_object.cladalSeqTotals = json.dumps([str(a) for a in cladalSeqAbundanceCounter])
            data_set_sample_object.save()


            dssList = []
            print('Populating sample {} with clade {} sequences'.format(samples[i], clade))
            if sum(countArray[i]) > 200:
                newCC = clade_collection(clade=clade, dataSetSampleFrom=data_set_sample_object)
                newCC.save()

            refSeqAbundanceCounter = defaultdict(int)
            for j in range(len(
                    nodes)):  # Only process this sample for this clade if there are enough seqs to make a cladeCollection
                abundance = countArray[i][j]
                if abundance > 0:
                    # I want to address a problem we are having here. Now that we have thorough checks to
                    # associate very similar sequences with indels by the primers to the same reference seq
                    # it means that multiple sequences within the same sample can have the same referenceseqs
                    # Due to the fact that we will in effect use the sequence of the reference seq rather
                    # than the dsss seq, we should consolidate all dsss seqs with the same reference seq
                    # so... we will create a counter that will keep track of the cumulative abundance associated with each reference_sequence
                    # and then create a dsss for each refSeq
                    refSeqAbundanceCounter[
                        reference_sequence.objects.get(id=nodeToRefDict[nodes[j]])] += abundance
            # > 200 associate a CC to the data_set_sample, else, don't
            # irrespective, associate a data_set_sample
            if sum(countArray[i]) > 200:
                for refSeq in refSeqAbundanceCounter.keys():
                    dss = data_set_sample_sequence(referenceSequenceOf=refSeq,
                                                   cladeCollectionTwoFoundIn=newCC,
                                                   abundance=refSeqAbundanceCounter[refSeq],
                                                   data_set_sample_from=data_set_sample_object)
                    dssList.append(dss)
                # Save all of the newly created dss
                data_set_sample_sequence.objects.bulk_create(dssList)
                # Get the ids of each of the dss and add create a string of them and store it as cc.footPrint
                # This way we can quickly get the footprint of the CC.
                # Sadly we can't get eh IDs from the list so we will need to re-query
                # Instead we add the ID of each refseq in the refSeqAbundanceCounter.keys() list
                newCC.footPrint = ','.join([str(refSeq.id) for refSeq in refSeqAbundanceCounter.keys()])
                newCC.save()
            else:
                for refSeq in refSeqAbundanceCounter.keys():
                    dss = data_set_sample_sequence(referenceSequenceOf=refSeq,
                                                   abundance=refSeqAbundanceCounter[refSeq],
                                                   data_set_sample_from=data_set_sample_object)
                    dssList.append(dss)
                # Save all of the newly created dss
                data_set_sample_sequence.objects.bulk_create(dssList)




    return



def main(pathToInputFile, dSID, numProc, screen_sub_evalue=False,
         full_path_to_nt_database_directory='/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload'):

    ### I will hard code in the mothur directory for the

    ##### CLEAN UP tempData FOLDER #### # this will need to be commented out if we are going to be starting from half way through an analysis, e.g. when debugging the MED work. Else the tempdata folder will be empty.
    # Delete the tempDataFolder and contents
    shutil.rmtree(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/tempData')))
    # recreate the tempDataFolder
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/tempData')),
                exist_ok=True)

    ############### UNZIP FILE, CREATE LIST OF SAMPLES AND WRITE stability.files FILE ##################

    dataSubmissionInQ = data_set.objects.get(id=dSID)
    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    if dataSubmissionInQ.initialDataProcessed == False:
        # create directory in which to write process files
        # wkd = r'{0}/{1}'.format(os.path.dirname(pathToInputFile), dSID)
        wkd = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/tempData/{}'.format(dSID)))

        try:
            os.makedirs(wkd)
        except FileExistsError:
            pass

        # Check to see if the files are already decompressed
        # If so then simply copy the files over to the destination folder
        compressed = True

        for file in os.listdir(pathToInputFile):

            print(file)
            if 'fastq.gz' in file or 'fq.gz' in file:
                # Then there is a fastq.gz already uncompressed in this folder
                # In this case we will assume that the seq data is not compressed into a mast .zip or .gz
                # Copy to the wkd
                compressed = False
                os.chdir('{}'.format(pathToInputFile))
                print(os.getcwd())
                # * asterix are only expanded in the shell and so don't work through subprocess
                # need to use the glob library instead
                 #https://stackoverflow.com/questions/13875978/python-subprocess-popen-why-does-ls-txt-not-work

                if 'fastq.gz' in file:
                    completedProcess = subprocess.run(['cp'] + glob.glob('*.fastq.gz') + [wkd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                elif 'fq.gz' in file:
                    completedProcess = subprocess.run(['cp'] + glob.glob('*.fq.gz') + [wkd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                break
        # Decompress the file to destination
        if compressed:
            extComponents = pathToInputFile.split('.')
            if extComponents[-1] == 'zip':  # .zip
                completedProcess = subprocess.run(["unzip", pathToInputFile, '-d', wkd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            elif extComponents[-2] == 'tar' and extComponents[-1] == 'gz':  # .tar.gz
                completedProcess = subprocess.run(["tar", "-xf", pathToInputFile, "-C", wkd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            elif extComponents[-1] == 'gz':  # .gz
                completedProcess = subprocess.run(["gunzip", "-c", pathToInputFile, ">", wkd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)


        # TODO see if we can make this simpler
        # Get rid of any dashes in the .gz filenames
        # Mothur does not allow '-' character in group names
        # When we use the filename anywhere but Mothur we can replace the tilde with the dash
        # Get string upto first underscore and consider it name
        # Check to see that it is not shared by more than two samples
        # Do this for each name. If it is shared then use the next underscore

        list_of_gz_files_in_wkd = [a for a in os.listdir(wkd) if '.gz' in a]

        # I think the simplest way to get sample names is actually to find what parts are common between all samples
        # well actually 50% of the samples so that we also remove the R1 and R2 parts.

        i = 1
        while 1:
            list_of_endings = []
            for file in list_of_gz_files_in_wkd:
                list_of_endings.append(file[-i:])
            if len(set(list_of_endings)) > 2:
                break
            else:
                i += 1
                # then this is one i too many and our magic i was i-1
        end_index = i-1

        list_of_names_non_unique = []
        for file in list_of_gz_files_in_wkd:
            list_of_names_non_unique.append(file[:-end_index])
        list_of_names = list(set(list_of_names_non_unique))

        if len(list_of_names) != len(list_of_gz_files_in_wkd)/2:
            print('Error in sample name extraction')
            return

        apples = 'asdf'


        # Make a batch file for mothur, set input and output dir and create a .file file

        mBatchFile = [
            r'set.dir(input={0})'.format(wkd),
            r'set.dir(output={0})'.format(wkd),
            r'make.file(inputdir={0}, type=gz, numcols=3)'.format(wkd)
        ]
        writeListToDestination(r'{0}/mBatchFile_makeFile'.format(wkd), mBatchFile)
        completedProcess = subprocess.run(['mothur', r'{0}/mBatchFile_makeFile'.format(wkd)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Convert the group names in the stability.files so that the dashes are converted to '[ds]',
        # So for the mothur we have '_'s. But for all else we convert these '_'s to dashes
        sampleFastQPairs = readDefinedFileToList(r'{0}/stability.files'.format(wkd))

        newstabilityFile = []
        for pair in sampleFastQPairs:
            pairComponenets = pair.split('\t')
            # I am going to use '[dS]' as a place holder for a dash
            # If I remember correctly each line of the stability file is a three column format with the first
            # column being the sample name. The second and third are the full paths of the .fastq.gz files

            newstabilityFile.append('{}\t{}\t{}'.format(
                pairComponenets[1].split('/')[-1][:-end_index].replace('-', '[dS]'),
                pairComponenets[1], pairComponenets[2]))



        writeListToDestination(r'{0}/stability.files'.format(wkd), newstabilityFile)
        sampleFastQPairs = newstabilityFile


        dataSubmissionInQ.workingDirectory = wkd
        dataSubmissionInQ.save()


        # Create data_set_sample instances
        listOfSamples = []
        for sampleName in list_of_names:
            print('Creating data_set_sample {}'.format(sampleName))
            # The cladalSeqTotals property of the data_set_sample object keeps track of the seq totals for the
            # sample divided by clade. This is used in the output to keep track of sequences that are not
            # included in cladeCollections
            emptyCladalSeqTotals = json.dumps([0 for cl in cladeList])
            dss = data_set_sample(name=sampleName, dataSubmissionFrom=dataSubmissionInQ, cladalSeqTotals=emptyCladalSeqTotals)  # We can do this using bulk_create
            listOfSamples.append(dss)
            # http://stackoverflow.com/questions/18383471/django-bulk-create-function-example
        smpls = data_set_sample.objects.bulk_create(listOfSamples)


    ################### CREATE CONTIGS SAMPLE BY SAMPLE #################


        print('Processing fastQ files')

        sampleFastQPairs = readDefinedFileToList(r'{0}/stability.files'.format(wkd))

        # Create the queues that will hold the sample information
        taskQueue = Queue()

        # Queue for output of successful and failed sequences
        outputQueue = Queue()

        # This will be a dictionary that we use to keep track of sequences that are found as matches when we do the
        # blast search against the symClade.fa database but that fall below the e value cutoff which is currently set
        # at e^-100. It will be a dictionary of sequence to number of samples in which the sequence was found in
        # the logic being that if we find sequences in multiple samples then they are probably genuine sequences
        # and they should therefore be checked against the full blast database to see if they match Symbiodinium
        # if they do then they should be put into the database.
        e_value_manager = Manager()
        e_value_multiP_dict = e_value_manager.dict()

        for contigPair in sampleFastQPairs:
            taskQueue.put(contigPair)


        for n in range(numProc):
            taskQueue.put('STOP')

        allProcesses = []
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()
        for n in range(numProc):
            p = Process(target=worker, args=(taskQueue, outputQueue, wkd, dSID, e_value_multiP_dict, dataSubmissionInQ.reference_fasta_database_used))
            allProcesses.append(p)
            p.start()

        for p in allProcesses:
            p.join()

        failedList = []

        outputQueue.put('STOP')

        for i in iter(outputQueue.get, 'STOP'):
            failedList.append(i)


        print('{0} out of {1} samples successfully passed QC.\n'
              '{2} samples produced erorrs'.format((len(sampleFastQPairs)-len(failedList)), len(sampleFastQPairs), len(failedList)))
        for contigPair in sampleFastQPairs:
            sampleName = contigPair.split('\t')[0].replace('[dS]', '-')
            if sampleName not in failedList:
                print('Sample {} processed successfuly'.format(sampleName))
        for sampleName in failedList:
            print('Sample {} : ERROR in sequencing reads. Unable to process'.format(sampleName))
        ### Here we should have all of the cladally sorted sequences printed out per sample.
        ### We now need to come up with a way of grouping them for MED.
        ### We also need to set initialDataProcessed to True
        dataSubmissionInQ.initialDataProcessed = True

        dataSubmissionInQ.save()

    # This function now performs the MEDs sample by sample, clade by clade.
    # The list of outputed paths lead to the MED directories where node info etc can be found
    MEDDirs = deuniqueFiles(dataSubmissionInQ.workingDirectory, dataSubmissionInQ.id, numProc)


    processMEDDataDirectCCDefinition_new_dss_structure(dataSubmissionInQ.workingDirectory, dataSubmissionInQ.id, MEDDirs)

    ###UNCOMMENT
    # dataSubmissionInQ.dataProcessed = True
    # dataSubmissionInQ.whenCompleted = datetime.now()
    dataSubmissionInQ.currentlyBeingProcessed = False
    dataSubmissionInQ.save()

    ### WRITE OUT REPORT OF HOW MANY SAMPLES WERE SUCCESSFULLY PROCESSED
    sampleList = data_set_sample.objects.filter(dataSubmissionFrom=dataSubmissionInQ)
    failedList = []
    for sample in sampleList:
        if sample.errorInProcessing:
            failedList.append(sample.name)

    readMeList = []
    sumMessage = '{0} out of {1} samples successfully passed QC.\n' \
                 '{2} samples produced erorrs'.format((len(sampleList) - len(failedList)), len(sampleList),
                                               len(failedList))
    print(sumMessage)
    readMeList.append(sumMessage)
    for sample in sampleList:

        if sample.name not in failedList:
            print('Sample {} processed successfuly'.format(sample.name))
            readMeList.append('Sample {} processed successfuly'.format(sample.name))
        else:
            print('Sample {} : {}'.format(sample.name, sample.errorReason))
    for sampleName in failedList:

        readMeList.append('Sample {} : ERROR in sequencing reads. Unable to process'.format(sampleName))
    writeListToDestination(pathToInputFile + '/readMe.txt', readMeList)

    # Here I also want to by default output a sequence drop that is a drop of the named sequences and their associated
    # sequences so that we mantain a link of the sequences to the names of the sequences
    sequence_drop_file = perform_sequence_drop()
    writeListToDestination(os.path.dirname(__file__) + '/dbBackUp/seq_dumps/seq_dump_' + str(datetime.now()), sequence_drop_file)



    ###### Assess and print out a fasta of the sequences that were found in multiple samples but were
    # below the e_value cutOff. Let's print these off in the directory that contains the data_set's
    # fastq.gz files.

    # At this point we know that they were found in multiple samples of the data_set which tells us that
    # they are very unlikely to be sequencing artefacts. But, we still need to check that they are Symbiodinium
    # To check that they are Symbiodinium we should do a blast and return something like 10 values and check to see
    # that 'Symbiodinium' is either the organism or in the name.

    # make fasta from the dict

    below_e_cutoff_dict = dict(e_value_multiP_dict)
    temp_count = 0
    fasta_out = []
    for key, value in below_e_cutoff_dict.items():
        if value > 2:
            # then this is a sequences that was found in three or more samples
            fasta_out.extend(['>sub_e_seq_{}_{}_{}'.format(dSID, temp_count, value), key])
            temp_count += 1
    writeListToDestination(wkd + '/blastInputFasta.fa', fasta_out)

    # we need to know what clade each of the sequences are
    # fastest way to do this is likely to run another blast on the symbiodinium clade reference dict
    ncbircFile = []
    db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))
    ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
    # write the .ncbirc file that gives the location of the db
    writeListToDestination("{0}/.ncbirc".format(wkd), ncbircFile)

    blastOutputPath = r'{}/blast.out'.format(wkd)
    outputFmt = "6 qseqid sseqid staxids evalue"
    inputPath = r'{}/blastInputFasta.fa'.format(wkd)
    os.chdir(wkd)

    # Run local blast
    # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
    completedProcess = subprocess.run(
        ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', dataSubmissionInQ.reference_fasta_database_used,
         '-max_target_seqs', '1', '-num_threads', '1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Read in blast output
    blast_output_file = readDefinedFileToList(r'{}/blast.out'.format(wkd))

    # now create a below_e_cutoff_seq to clade dictionary
    sub_e_seq_to_clade_dict = {a.split('\t')[0]: a.split('\t')[1][-1] for a in blast_output_file}

    # print out the fasta with the clades appended to the end
    fasta_out_with_clade = []
    for line in fasta_out:
        if line[0] == '>':
            fasta_out_with_clade.append(line + '_clade' + sub_e_seq_to_clade_dict[line[1:]])
        else:
            fasta_out_with_clade.append(line)

    # this will return a new fasta containing only the sequences that were 'Symbiodinium' matches
    # we can then output this dictionary
    writeListToDestination(pathToInputFile + '/below_e_cutoff_seqs_{}.fasta'.format(dSID), fasta_out_with_clade)

    ##### CLEAN UP tempData FOLDER #####
    # Delete the tempDataFolder and contents
    shutil.rmtree(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/tempData')))
    # recreate the tempDataFolder
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'SymPortal_Data/tempData')),
                exist_ok=True)

    # write out whether there were below e value sequences outputted.
    print('WARNING: {} sub_e_value cut-off sequences were output'.format(len(fasta_out_with_clade)/2))
    if screen_sub_evalue:
        print('These will now be automatically screened to see if they contain Symbiodinium sequences.')
        print('Screening sub e value sequences...')

        apples = 'pears'

        symportal_framework_object = symportal_framework.objects.get(id=1)
        preious_reference_fasta_name = symportal_framework_object.latest_reference_fasta
        required_sample_support = symportal_framework_object.required_sub_e_value_seq_support_samples
        required_symbiodinium_blast_matches = symportal_framework_object.required_sub_e_value_seq_support_blast_symbiodinium
        next_reference_fasta_iteration_id = symportal_framework_object.next_reference_fasta_iteration
        new_reference_fasta_name, num_additional_sequences, new_ref_fasta_location = screen_sub_e_value_sequences(dSID, pathToInputFile, iteration_id=next_reference_fasta_iteration_id, seq_sample_support_cut_off=required_sample_support, previous_reference_fasta_name=preious_reference_fasta_name, required_symbiodinium_matches=required_symbiodinium_blast_matches, full_path_to_nt_database_directory=full_path_to_nt_database_directory)

        print('Done\n')

        if new_reference_fasta_name:

            print('WARNING: {} Symbiodinium sequences were found in those discarded due to e value cutoffs.'.format(num_additional_sequences))
            print('A new reference fasta has been created that contains these new sequences as well as those that '
                  'were contained in the previous version of the reference fasta.')
            print('This new reference fasta is called: {}'.format(new_reference_fasta_name))
            print('It has been output to the following location: {}'.format(new_ref_fasta_location))
        else:
            print('None of the e value discarded sequences returned matches for Symbiodinium when run against '
                  'the nt database.\nHappy days!')
    else:
        print('To screen these sequences for possible symbiodinium sequences please set screen_sub_evalue '
              'to True and provide a directory that contains the NCBI nt database')

def screen_sub_e_value_sequences(ds_id, data_sub_data_dir, iteration_id, seq_sample_support_cut_off, previous_reference_fasta_name, required_symbiodinium_matches, full_path_to_nt_database_directory):
    # we need to make sure that we are looking at matches that cover > 95%
    # this is probably the most important point. We can then decide what percentage coverage we want
    # perhaps something like 60%.
    # we then need to see if there is a 'Symbiodinium' sequence that matches the query and all of these
    # requirements. If so then we consider the sequence to be Symbiodinium
    # TODO make sure that we have metrics that show how many sequences were kicked out for each iterarion that we
    # do the database update.
    # We should write out the new database with an iteration indicator so that we can keep track of the progress of the
    # database creations. We can then run the database submissions using specific iterations of the symclade dataase
    # we can name the data_set that we do so that they can link in with which database iteration they are using

    # we can work with only seuqences that are found above a certain level of support. We can use the
    # seq_sample_support_cut_off for this.


    # Write out the hidden file that points to the ncbi database directory.
    ncbircFile = []
    # db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))

    db_path = full_path_to_nt_database_directory
    ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
    writeListToDestination("{}.ncbirc".format(data_sub_data_dir + '/'), ncbircFile)


    # Read in the fasta files of below e values that were kicked out.
    fasta_file = readDefinedFileToList('{}/below_e_cutoff_seqs_{}.fasta'.format(data_sub_data_dir, ds_id))
    fasta_file_dict = createDictFromFasta(fasta_file)

    # screen the input fasta for sample support according to seq_sample_support_cut_off
    screened_fasta = []
    for i in range(len(fasta_file)):
        if fasta_file[i][0] == '>':
            if int(fasta_file[i].split('_')[4]) >= seq_sample_support_cut_off:
                screened_fasta.extend([fasta_file[i], fasta_file[i + 1]])

    # write out the screened fasta so that it can be read in to the blast
    # make sure to reference the sequence support and the iteration
    path_to_screened_fasta = '{}/{}_{}_{}.fasta'.format(data_sub_data_dir,'below_e_cutoff_seqs_{}.screened'.format(ds_id), iteration_id, seq_sample_support_cut_off)
    screened_fasta_dict = createDictFromFasta(screened_fasta)
    writeListToDestination(path_to_screened_fasta, screened_fasta)

    # Set up environment for running local blast
    blastOutputPath = r'{}/blast_{}_{}.out'.format(data_sub_data_dir, iteration_id, seq_sample_support_cut_off)
    outputFmt = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"
    # inputPath = r'{}/below_e_cutoff_seqs.fasta'.format(data_sub_data_dir)
    os.chdir(data_sub_data_dir)

    # Run local blast
    # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
    completedProcess = subprocess.run(
        ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', path_to_screened_fasta, '-db', 'nt',
         '-max_target_seqs', '10', '-num_threads', '20'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Read in blast output
    blast_output_file = readDefinedFileToList(r'{}/blast_{}_{}.out'.format(data_sub_data_dir, iteration_id, seq_sample_support_cut_off))

    # blastDict = {a.split('\t')[0]: a.split('\t')[1][-1] for a in blastOutputFile}
    apples = 'asdf'
    # Now go through each of the results and look to see if there is a result that has > 95% coverage and has >60%
    # match and has symbiodinium in the name.
    # if you find one then add the name of this seq to the reference db

    # create a dict that is the query name key and a list of subject return value
    blast_output_dict = defaultdict(list)
    for line in blast_output_file:
        blast_output_dict[line.split('\t')[0]].append('\t'.join(line.split('\t')[1:]))

    verified_sequence_list = []
    for k, v in blast_output_dict.items():
        sym_count = 0
        for result_str in v:
            if 'Symbiodinium' in result_str:
                percentage_coverage = float(result_str.split('\t')[4])
                percentage_identity_match = float(result_str.split('\t')[3])
                if percentage_coverage > 95 and percentage_identity_match > 60:
                    sym_count += 1
                    if sym_count == required_symbiodinium_matches:
                        verified_sequence_list.append(k)
                        break

    # We only need to proceed from here to make a new database if we have sequences that ahve been verified as
    # Symbiodinium
    if verified_sequence_list:
        # here we have a list of the Symbiodinium sequences that we can add to the reference db fasta
        new_fasta = []
        for seq_to_add in verified_sequence_list:
            new_fasta.extend(['>{}'.format(seq_to_add), '{}'.format(screened_fasta_dict[seq_to_add])])

        # now add the current sequences
        previous_reference_fasta = readDefinedFileToList('{}/{}'.format(os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')), previous_reference_fasta_name))
        # we need to check that none of the new sequence names are found in
        new_fasta += previous_reference_fasta

        # now that the reference db fasta has had the new sequences added to it.
        # write out to the db to the database directory of SymPortal
        full_path_to_new_ref_fasta_iteration = '{}/symClade_{}_{}.fa'.format(os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB')), iteration_id, seq_sample_support_cut_off)
        writeListToDestination(full_path_to_new_ref_fasta_iteration, new_fasta)

        # now update the SymPortal framework object
        symportal_framework_object = symportal_framework.objects.get(id=1)
        symportal_framework_object.latest_reference_fasta = 'symClade_{}_{}.fa'.format(iteration_id, seq_sample_support_cut_off)
        symportal_framework_object.next_reference_fasta_iteration += 1
        symportal_framework_object.save()

        # run makeblastdb
        completed_process = subprocess.run(['makeblastdb','-in', full_path_to_new_ref_fasta_iteration,'-dbtype' ,'nucl', '-title', 'symClade_{}_{}'.format(iteration_id, seq_sample_support_cut_off)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        return 'symClade_{}_{}.fa'.format(iteration_id, seq_sample_support_cut_off), len(verified_sequence_list), full_path_to_new_ref_fasta_iteration
    else:
        return False, 0, False

def perform_sequence_drop():
    # this will simply produce a list
    # in the list each item will be a line of text that will be a refseq name, clade and sequence
    output_list = []
    for ref_seq in reference_sequence.objects.filter(hasName=True):
        output_list.append('{}\t{}\t{}'.format(ref_seq.name, ref_seq.clade, ref_seq.sequence))
    return output_list