console.log('snpApp :-)')

// Assemble application under div id=snpApp

// 1. Get SNP textArea going
snpApp={

	disp:function(x){
		//console.log(x);
		var li = document.createElement('li');
		li.innerHTML=x; // note innerHTML allows HTML tagging, not as plain as textContent
		jmat.gId('disp').appendChild(li);
		return li;
	},

	start:function(){
		// add text area to collect SNP string
		var calc = jmat.gId('snpCalc')
		calc.onclick=function(evt){
			jmat.gId('disp').innerHTML='<hr>';
			snpApp.disp('starting ...');
			snpApp.parseSnp(evt);
			snpApp.disp('SNPs: '+snpApp.dt.snps[0].length);
			var n=snpApp.dt.rows.length;
			snpApp.disp('Strains:');
			for (i=0;i<n;i++){snpApp.disp(i+'/'+n+' - '+snpApp.dt.rows[i]+' ['+snpApp.dt.y[i]+']')}
			// identify partition table
			snpApp.disp('computing partition table ...')
			snpApp.dt = snpCalc.partitions(snpApp.dt);
			snpApp.disp('ranked partitions, SNPs matching <input type=text id="SNPquery">');
			I = jmat.array2mat(snpApp.dt.partition); // seprating binary indexes from ranksum values
			J=jmat.sort(I[1].map(function(ii){return -ii})); // sorting ranksums, descending
			snpApp.dt.Ind=J[1];
			snpApp.dt.partitionCodeList=I[0];
			snpApp.dt.rankSumList=I[1];
			snpCalc.assignSnps();//Assigns SNPs to partitions
			var m = snpApp.dt.Ind.length;
			for(j=0;j<m;j++){ // look carefully here, it explain how to index ranked results
				var ii = snpApp.dt.Ind[j];
				li = snpApp.disp(j+'/'+m+' : '+snpApp.dt.partitionCodeList[ii]+' --> '+snpApp.dt.rankSumList[ii]+' ('+snpApp.dt.matchedSNPs[ii].length+')')
				li.id = snpApp.dt.partitionCodeList[ii];
				li.onclick=function(x){snpApp.showSnps(x)};
			}

		}
		//snpApp.load(); // this should be activated by Select
	},
	
	showSnps:function(x){ // mouse click on closed partition
		var li = x.target;
		var i = snpApp.dt.partition[li.id]; //partition being opened
		var jj=snpApp.dt.matchedSNPs[i]; // SNPs that match it
		jj.map(function(j,ij){
			// retrieve filter
			var q=jmat.gId('SNPquery').value;
			if (q.length==0){q='.'}
			q=RegExp(q,'g');
			if (snpApp.dt.cols[j].match(q)){
				var i = snpApp.dt.partition[li.id]; //partition being opened
				var snp = snpApp.dt.snps.map(function(s){return s[j]});
				snp=jmat.array2str(snp,'0').replace(/0/g,'');
				var p =document.createElement('p');
				p.className=li.id;
				p.innerHTML=(ij+1)+'. ['+snpApp.dt.matchedVars[i][ij]+'] <span style="color:red">'+snpApp.dispSnp(snp)+'</span> '+snpApp.dt.cols[j];
				p.style.color='green';
				li.appendChild(p);
			}
		});
		li.onclick=function(x){snpApp.hideSnps(x)};
	},
	
	dispSnp:function(x){// displaying the snp string
		return x;
	},
	
	hideSnps:function(x){ // mouse click on opened partition
		var pp = document.getElementsByClassName(x.target.id);
		while(pp.length>0){x.target.removeChild(pp[0])}
		x.target.onclick=function(x){snpApp.showSnps(x)};
	},
	
	loadSelect:function(evt){
		var dt={}, fname=evt.target.value;
		snpApp.disp('Loading '+fname);
		jmat.load(fname);
	},

	load:function(evt){ // load existing SNP set
		var dt={};
		jmat.load('dtPneumopath.js');		
	},

	dtCallback:function(dt){ // included in loaded file with dt
		//var org = jmat.gId('organism');
		//org.value=jmat.array2str(dt.rows);
		//org.rows=dt.rows.length
		var snp = jmat.gId('SNPs');
		var n=dt.rows.length;
		snp.rows=n;
		var snps=[];
		for (var i=0;i<n;i++){
			snps[i]=dt.rows[i]+'\t'+dt.snps[i];
		}
		snp.value=jmat.array2str(snps);
		var gen = jmat.gId('genotypes');
		gen.value=jmat.array2str(dt.cols);
		gen.rows=n;
		//console.log(dt);
		snpApp.parseSnp()
		snpApp.parseGenotypes()
	},

	parseSnp:function(evt){
		var linhas = jmat.gId('SNPs').value.split('\n');
		snpApp.disp('SNP parsing ...');
		if(!snpApp.dt){
			snpApp.dt={snps:[],rows:[],cols:[],y:[]};
		}
		var n = linhas.length;
		for (var i=0;i<n;i++){
			var linha = linhas[i].split('\t');
			snpApp.dt.snps[i]=linha[1];
			snpApp.dt.rows[i]=linha[0];	
		}
		jmat.gId('organisms').innerHTML=jmat.array2str(snpApp.dt.rows,'<br>');
		jmat.gId('phenotypes').rows=n;
		
	},
	
	parseGenotypes:function(evt){
		if(!snpApp.dt){snpApp.dt={snps:[],rows:[],cols:[],y:[]};}
		var linhas = jmat.gId('genotypes').value.split('\n');
		snpApp.disp('Genotype parsing ...');
		for (var i=0;i<linhas.length;i++){snpApp.dt.cols[i]=linhas[i];}
	}
};

snpCalc={
	partitions:function(dt){ // identifies p value for all possible partition profiles
		var linhas = jmat.gId('phenotypes').value.split('\n');
		snpApp.disp('Phenotype parsing ...');
		var j=0;
		for (var i=0;i<linhas.length;i++){
			snpApp.dt.y[j]=linhas[j].split(',').map(function(x){return JSON.parse(x)});
			j=j+1;	
		}
		var n=dt.rows.length, pt=[];
		for(var i=1;i<=Math.pow(2,n)-2;i++){
			pt[i-1]=jmat.dec2bin(i,n);
		}
		dt.partition=[];
		m=pt.length;
		for(var i=0;i<m;i++){
			var x = [0]; // just to seed the concatenation
			var y = [0];
			for(var j=0;j<n;j++){
				if(pt[i][j]=='1'){x=jmat.cat(x,dt.y[j])}
				else{y=jmat.cat(y,dt.y[j])}
			}
			dt.partition[pt[i]]=jmat.ranksum(x.slice(1),y.slice(1)); //note seed sliced off
		}
		return dt
	},
	assignSnps:function(){ // assign SNPs to the first n partitions found
		// create an array of SNPs per position
		//var snps=dt.snps.map(function(s){
		//	var si='';
		//	for (var i=0;i<=n;i++){
		//		si+=s[dt.Ind[i]];				
		//	}	
		//	return si;				
		//})
		var j=0;for(i in snpApp.dt.partition){snpApp.dt.partition[i]=j;j++} // number partitions
		var snps = jmat.transpose(snpApp.dt.snps.map(function(x){return x.split('')}));
		snpApp.dt.matchedSNPs=snpApp.dt.partitionCodeList.map(function(x){return []}); // creates array of the SNPs matched with each partition
		snpApp.dt.matchedVars=snpApp.dt.partitionCodeList.map(function(x){return []}); // the variants of each match
		snps.map(function(s,ii){ // ith SNP
			//var j=dt.Ind[jj];
			var u = jmat.unique(s);
			u.map(function(ui){
				//determine the partition string
				var pi=jmat.array2str(s.map(function(si){return ''+(si==ui)*1}),'.').replace(/\./g,'');
				var j=snpApp.dt.partition[pi]; 
				//console.log(j);
				var k=snpApp.dt.matchedSNPs[j].length;
				snpApp.dt.matchedSNPs[j][k]=ii; // partition #j is matched by SNP #ii
				snpApp.dt.matchedVars[j][k]=ui; // by variant ui
				//var pr = document.createElement('p');
				//pr.innerHTML=ui+': ['+jmat.array2str(s,'.').replace(/\./g,'')+'] SNP#'+j+' '+snpApp.dt.cols[j];
				//pr.className='part'+pi+' snp'+j;
				//jmat.gId(pi).appendChild(pr);
			})
		})
	}
}

snpApp.start();

