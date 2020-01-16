#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class Pipeline:
    
    def __init__(self,organism,gene,prot,email="paulo-carvalhais@hotmail.com",threshold=0.001):
        self._email=email
        self.lookup_dic={}
        self.threshold=threshold
        Entrez.email=self._email
        self.organism=organism
        self.gene=gene
        self.prot_id=prot
        self.gene_id=''
    
    def extract_abstracts(self): ## Extrair Abstracts de artigos no pubmed
        '''
        Função que procura um Gene e um Organismo contra a base de dados do pubmed.
        Cria um ficheiro de texto para cada artigo encontrado, com o seu Id, Titlo e Abstract
        para fácil análise.
        '''
        handle = Entrez.esearch(db="pubmed", term=self.organism+"[Orgn] AND "+self.gene+"[Gene]",retmode='text') ##Base de dados a procurar, especificando o organismo e o Gene
        record= Entrez.read(handle) 
        id_lists=record['IdList'] ## Colecção de Ids extraídos da procura
        handle.close()
    
        for i in id_lists:
            handle = Entrez.efetch(db='pubmed',id=i, rettype='Medline',retmode='text')
            record= Medline.read(handle)
            file=open(self.gene+'_'+i+'.txt','w')
    
            file.write('Pubmed_Id: '+ i+'\n')
            file.write('Title: \n'+ record['TI']+'\n')
            file.write('Abstract: \n'+ record['AB']+'\n')
            
            file.close()
            
    def get_gene_id(self):
        '''
        Função que procura o id do Gene no organismo especificados
        Grava um ficheiro no sistema com informação em formato xml
        '''
        handle = Entrez.esearch(db="gene", term=self.organism+"[Orgn] AND "+self.gene+"[Gene]",retmode='text') ##Base de dados a procurar, especificando o organismo e o Gene
        record= Entrez.read(handle)
        gi_list=record['IdList']
        handle=Entrez.efetch(db='gene',id=gi_list[0],retmode='xml')
        with open(self.gene+'.txt','w') as out_handle:
            out_handle.write(handle.read())
        self.gene_id=gi_list[0]
        return self.gene_id
        
    def get_genbank_records(self):
        '''
        Função que procura ficheiros de formato GenBank relacionados com o organismo e o Gene desejados.
        Cria vários ficheiros GenBank
        '''
        handle = Entrez.esearch(db='nucleotide', term=self.organism+"[Orgn] AND "+self.gene+"[Gene]",usehistory='y') ##Pesquisa especifica na base de dados
        record = Entrez.read(handle)
        webenv = record["WebEnv"]
        query_key = record["QueryKey"]
        gi_list=record['IdList'] ## Retirar lista de Ids dos ficheiros GenBank
        handle.close()
        dic={}
        records=[]
        genomic_ids=[]
        for id_ in gi_list:
            handle=Entrez.efetch(db='nucleotide',id=id_,rettype='gb',retmode='text',webenv=webenv,query_key=query_key)
            record=SeqIO.read(handle,'gb')
            records.append(record)
            handle.close()
            for rec in records:
                for feat in rec.features:
                    for key, val in feat.qualifiers.items():
                        if 'genomic DNA' in val:
                            genomic_ids.append(rec.id)
                        
            for rec in records:
                if rec.id in genomic_ids:
                    file=open(self.gene+'_GBfile_'+id_+'.txt','w')
                    file.write(Entrez.efetch(db='nucleotide',id=id_,rettype='gb',retmode='text').read())
                    file.close
                for feat in rec.features:
                    if feat.type=='CDS':    
                        dic[rec.name]=(feat.qualifiers['protein_id'],feat.qualifiers['translation'])

              
    def blast(self,program='blastn',base='nt',evalue=0.001,align=50,matrix='BLOSUM62'):
        '''
        Função que corre o blast de sequencias de dna contra a base de dados pretendida, tendo como parametros default
        o programa-blastn
        a matrix- BLOSUM62
        o e-value=0.001
        o nº de alinhamentos maximos=100
        '''
        
        result_handle=NCBIWWW.qblast(program,'nt',self.gene) ### Bio.SearchIO
        with open(self.gene+'_blast.xml','w') as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        
        hit_id=[]
        result_handle = open(self.prot_id+'_blast.xml')
        blast_records = NCBIXML.parse(result_handle)
        f = open(self.prot_id+'_seqs.fasta', 'w+')
        e_value_thresh = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        hit_id.append(alignment.hit_id)
                        print("****Alignment****")
                        print("sequence:", alignment.title)
                        print("length:", alignment.length)
                        print("e value:", hsp.expect)
                        print(hsp.query[0:75] + "...")
                        print(hsp.match[0:75] + "...")
                        print(hsp.sbjct[0:75] + "...")
                        
        hit_id = list(set(hit_id))                  
        for ind in hit_id:               
            handle = Entrez.efetch(db ="protein", rettype ='fasta', retmode ="text", id=ind)
            seq_record = SeqIO.read(handle, 'fasta')
            handle.close()
            print('>'+seq_record.description, file=f)
            print(seq_record.seq, file=f)
            print('', file=f)

        file.close()
        
    def blast_p(self,progr='blastp',base='swissprot',evalue=0.001,align=50,matrix='BLOSUM62'):
        '''
        Função que corre o blast de sequencias de dna contra a base de dados pretendida, tendo como parametros default
        a base de dados-swissprot
        o programa-blastp
        a matrix- BLOSUM62
        o e-value=0.001
        o nº de alinhamentos maximos=50
        '''
        result_handle=NCBIWWW.qblast(program=progr,db=base,seq=self.prot_id,ncbi_gi=True,expect=evalue,alignments=align,matrix_name=matrix) ### Bio.SearchIO
        with open(self.prot_id+'_blast.xml','w') as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        
        hit_id=[]
        result_handle = open(self.prot_id+'_blast.xml')
        blast_records = NCBIXML.parse(result_handle)
        f = open(self.prot_id+'_seqs.fasta', 'w+')
        e_value_thresh = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        hit_id.append(alignment.hit_id)
                        print("****Alignment****")
                        print("sequence:", alignment.title)
                        print("length:", alignment.length)
                        print("e value:", hsp.expect)
                        print(hsp.query[0:75] + "...")
                        print(hsp.match[0:75] + "...")
                        print(hsp.sbjct[0:75] + "...")
                        
        hit_id = list(set(hit_id))                  
        for ind in hit_id:               
            handle = Entrez.efetch(db ="protein", rettype ='fasta', retmode ="text", id=ind)
            seq_record = SeqIO.read(handle, 'fasta')
            handle.close()
            print('>'+seq_record.description, file=f)
            print(seq_record.seq, file=f)
            print('', file=f)
        
    def multiple_alignment(self):
        '''
        Função que realiza o alinhamento multiplo com as melhores sequencias 
        do blast corrido anteriormente, e escolhidas pela função anterior 
        '''
        cmdline=ClustalwCommandline(r'C:\Program Files (x86)\ClustalW2\clustalw2',infile=self.prot_id+'_blast.fasta')
        cmdline()
        
    def create_tree(self):
        '''
        Função que cria a árvore filogenética com base no alinhamento múltiplo
        '''
        infile=(self.prot_id+'_seqs.dnd')
        tree = Phylo.read(infile, "newick")
        Phylo.draw_ascii(tree)
        return tree
        
    def create_tree_xml(self):
        '''
        Cria a arvore em formato phyloxml de modo a poder ser exportada para outros programas de desenho de arvores filogenéticas
        '''
        outfile=(self.prot_id+'_tree.xml')
        tree1=self.create_tree()
        tree_xml = tree1.as_phyloxml()            
        Phylo.write(tree_xml, outfile, 'phyloxml')
            
    def get_Uniprot_rec(self):
        '''
        Através do id da proteina codificada pelo gene, fornecido pelo utilizador, vai fazer o download das informações relativas
        à proteina encontradas na página correspondente na UniProt
        '''
        handle = ExPASy.get_sprot_raw(self.prot_id)
        url=handle.url
        url=url.replace('txt','xml')
        response = requests.get(url)
        with open('Uniprot'+self.prot_id+'.xml','wb') as file:
            file.write(response.content)
        
    def extract_info_uniprot(self):
        '''
        Através do ficheiro guardado anteriormente, analisa e extrai do mesmo informação que consideramos relevante em relação 
        à proteina e guarda-as num ficheiro txt
        '''
        lista_xrefs=[] # Lista de entradas de bases de dados referentes ao gene e ids correspondentes 
        PDB_tuples=[] ## Lista de tuplos da base de dados PDB e os ids que lhe correspondem
        PDB_ids=[] ### Lista de Ids da PDB
        info=SeqIO.read('Uniprot'+self.prot_id+'.xml','uniprot-xml')

        file=open('Uniprot_info_'+self.prot_id+'.txt','w')
        try:
            name=str(info.name)
            print('Nome:\n'+ name+'\n')
            file.write('Nome:\n'+ name+'\n')
        except:
            print('Ficheiro não contém informação de nome')
        try:
            alt_name=str(info.annotations['alternativeName_fullName'])
            print('Nome Alternativo:\n'+ alt_name+'\n')
            file.write('Nome Alternativo:\n'+ alt_name+'\n')
        except:
            print('Ficheiro não contém informação de nome alternativo')
        try:
            seq=str(info.seq)
            print('Sequencia:\n'+seq+'\n')
            file.write('Sequencia:\n'+seq+'\n')
        except:
            print('Ficheiro não contém informação de sequência ')
        try:
            function=str(info.annotations['comment_function'])
            print('Função:\n'+function+'\n')
            file.write('Função:\n'+function+'\n')
        except:
            print('Ficheiro não contém informação de função')
        try:
            sub_location=str(info.annotations['comment_subcellularlocation_location'])
            print('Localização:\n'+ sub_location+'\n')
            file.write('Localização:\n'+ sub_location+'\n')
        except:
            print('Ficheiro não contém informação de localização subcelular ')
        try:
            taxo=str(info.annotations['taxonomy'])
            print('Taxonomia:\n'+taxo+'\n')
            file.write('Taxonomia:\n'+taxo+'\n')
        except:
            print('Ficheiro não contém informação de taxonomia')  
        try:
            ptm=str(info.annotations['comment_PTM'])
            print('PTM:\n'+ptm+'\n')
            file.write('PTM:\n'+ptm+'\n')
        except:
            print('Ficheiro não contém informação de PTM')
        try:
            feats=info.features
            file.write('Features:\n')
            print('Features:\n')
            for i in feats:
                print(str(i)+'\n')
                file.write(str(i)+'\n')
        except:
            print('Ficheiro não contém informação de features')
        try:
            lista_xrefs=[] # Lista de entradas de bases de dados referentes ao gene e ids correspondentes 
            lista_xrefs.append(str(info.dbxrefs))
            print('CrossReferences:\n'+str(lista_xrefs)+'\n')
            file.write('CrossReferences:\n'+str(lista_xrefs)+'\n')
        except:
            print('Ficheiro não contém informação de cross-reference de bases de dados')
        file.close()

        for i in info.dbxrefs: ### Cria tuplos com base nas cross references de bases de dados
            a=tuple(i.split(':'))
            lista_xrefs.append(a)

        for i in lista_xrefs: ### Extrai apenas os tuplos com informação da PDB
            if i[0]=='PDB':
                PDB_tuples.append(i)

        for i in PDB_tuples: ### Extrai apenas os ids
            PDB_ids.append(i[1])
        PDB_error='A Proteína não possui ids de estrutura 3D na base de dados PDB'
        
        if len(PDB_ids)!=0:
            return ('Lista de IDs da PDB') + str(PDB_ids)
        else:
            return PDB_error
        
    def main(self):
        '''
        Corre todas as funções anteriores em sequência, criando uma pipeline automatizada de procura de informação
        '''
        self.extract_abstracts()
        self.get_gene_id()
        self.get_genbank_records()
        self.blast_p()
        self.multiple_alignment()
        self.create_tree_xml()
        self.get_Uniprot_rec()
        self.extract_info_uniprot()
    
if __name__ == '__main__':
    slc=Pipeline('Homo Sapiens','SLC35B2','Q8TB61')
    fox=Pipeline('Homo Sapiens','FOXC2','Q99958')
    zbt=Pipeline('Homo Sapiens','ZBTB7A','O95365')
    #slc.main()
    #fox.main()
    #zbt.main()

