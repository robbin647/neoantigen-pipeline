import Oviz from "crux";
import { editorConfig, editorRef } from "./editor";
import template from "./template.bvt";

import { ComplexBoxplot, processBoxData } from "oviz-components/complex-boxplot";
import { GridPlot } from "oviz-components/grid-plot";
import { EditText } from "oviz-components/edit-text";
import {register} from "page/visualizers";
import { rankDict, sortByRankKey } from "utils/bio-info";
import { registerEditorConfig } from "utils/editor";
import { generatePath, plotData } from "./data";
import * as d3 from "d3";
import * as ds from "d3-sankey";
import * as TextSize from "crux/dist/utils/text-size";
import { TableSimplePlugin } from "bootstrap-vue";
import { isFunctionOrConstructorTypeNode, isTemplateExpression, nodeModuleNameResolver } from "typescript";
import { tooltip } from "crux/dist/utils";
import { co2servermap,imageMap } from "../overview_species_tree/commonInfo"
import { Text } from "crux/src/element";
import { max } from "crux/dist/utils/math";
import { replace } from "utils/replace";
import axios from "axios";
// import { parse } from "utils/newick";

export class GeneGome extends Oviz.Component {

    public newdata: Node; 
    public newconverlist:{data:any[],column:any[]};

    public render() {
        return this.t`${template}`;
    }
    
    public willRender() {
    
        // // this.loadProteinData()
        // this.chosenGenesymbol != ""? this.genetree = this.childrentree:null
        // // // this.$v.size.height = 200
        // // this.chosenGenesymbol != ""? this.treeheight = 200:null

        // if(this.chosenGenesymbol != ""){
        //     // let data = 
        //     let data = this.orifamilyinfo
        //     console.log("gene tree:",data)
        //     let geneSymbol_list = Array.from(new Set(data.map(x=>x["gene_symbol"]))).filter(x=>x!= "").filter(x=>x!=undefined)
        //     let proteinid_list = Array.from(new Set(data.map(x=>x["protein_id"]))).filter(x=>x!= "").filter(x=>x!=undefined)
        //     let species_list = Array.from(new Set(data.map(x=>x["species"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "").filter(x=>x!=undefined)
        //     let maxsymlen = calcuSymbollen(geneSymbol_list)
        //     let infoData = {}
        //     data.forEach(element => {
        //         infoData[element.protein_id] = element 
        //     });

        //     let specieslist = Array.from(new Set(data.map(x=>x["species"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "")
        
        //     this.infoData = infoData
        //     console.log("this.infoData:",this.infoData)
        //     this.proteinid_list = proteinid_list
        //     this.specieslist = specieslist
        //     this.maxsymlen = maxsymlen
        //     this.geneSymbol_list = geneSymbol_list
        //     this.textSymbol_list = this.geneSymbol_list.map(x=>x.replace("gene_symbol:",""))
        //     this.chosensymlist = this.geneSymbol_list.slice(0,5)
        //     this.co2servermap = co2servermap

        // let childnodeslist = getList([this.genetree],this.infoData)
        // this.childnodeslist = childnodeslist

        // let colorMap = {}
        // species_list.forEach((ele,index) => {
        //     const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
        //     colorMap[ele]  = color.string
        // });
        // this.colorMap = colorMap

        // this.showproteincolorMap = false
        // this.proteincolorMap = {}
        // this.showspeciescolorMap =false
        // this.speciescolorMap = {}
        // this.showspeciescolorMap =false
        // this.partscolorMap = {}
        
        // this.showTitlelist = []
        // this.titlesep = getMaxlength(proteinid_list) + 10
        // this.geneMaxlen = getMaxlength(geneSymbol_list)


        // // this.chosenGenesymbol = ""
        // this.selectgeneSymbol_list = {select:this.geneSymbol_list}
        // //this.data.geneSymbol_list //所有的genelist

        // //yuan setup
        
        
        // this.geneBlockcolormap = {}
        // this.geneSymbol_list.forEach((element,index) => {
        //     const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
        //     this.geneBlockcolormap[element] = color.string
        // });


        // let genegroupData = groupBy(this.childnodeslist, (childnode) => {
        //         return childnode.genesymbol
        // })
        // let geneblockwidths = {}
        // Object.keys(genegroupData).forEach((item,index)=>{
        //     let temparr = []
        //     genegroupData[item].forEach(ele => {
        //         console.log("this.infoData[ele.name]:",this.infoData[ele.name])
        //         console.log("ele.name:",ele.name)
        //         let spiname = this.infoData[ele.name].species.replace("Lys_pep_need_to_confirm_","")
        //         let spitext = (co2servermap["server-sci"][spiname] + " / " +co2servermap["server-com"][spiname]).replace(/_/g," ") + "" + this.leafSeparator + this.infoData[ele.name].gene_symbol
        //         temparr.push(TextSize.measuredTextSize(spitext,10).width)
        //     });
        //     geneblockwidths[item.replace(/\"/g,"")] = max(temparr)
        // })
        // this.geneblockwidths = geneblockwidths

        // let protein2genelist = this.childnodeslist.map(x=>x["genesymbol"])
        // const findLongestSequence = (arr = []) => {
        //     const res = arr.reduce((acc,val,ind) => {
        //         if(acc.length && acc[acc.length-1][0] === val){
        //             acc[acc.length-1].push(val);
        //         }else{
        //             acc.push([val]);
        //         };
        //         return acc;
        //     },[])
        //     return res;
        // }
        // let geneblocks = findLongestSequence(protein2genelist)
        // this.geneblocks = geneblocks

        // this.treeheight = this.childnodeslist.length * 12 
        // }

    }
    
    private async loadProteinData() {
        const dataL = await axios.get(`/data/demo/genefamily_tree/childrenTree.tree.newick.txt`);
        this.newdata = parse(dataL)
        const convertlist = await axios.get(`/data/demo/genefamily_tree/childrenprotein2geneall.txt`);
        this.newconverlist = convetTSV(convertlist)

        
    }

}

export function convetTSV(oridata){
    let convertdata = {data:[],column:[]}
    let columname = oridata.data.split("\n")[0].split("\t")
    oridata.data.split("\n").forEach((element,index) => {
        let singlearray = {}
        element.split("\t").forEach((eachchr,chrindex) => {
            singlearray[columname[chrindex]] = eachchr
        });
        convertdata.data.push(singlearray)
    });
    convertdata.column = columname
    convertdata.data = convertdata.data.splice(1)

    return convertdata
}

let temp = []
export function getList(tree,infoData) {
    for (let item in tree) {
      if (tree[item].children) {
        getList(tree[item].children,infoData)
      } else {
        // console.log("infoData[tree[item].name].gene_symbol:",tree[item].name)
        infoData[tree[item].name]!= undefined ?
        tree[item]["genesymbol"] = infoData[tree[item].name].gene_symbol :null
          
        temp.push(tree[item])
      }
    }
    return temp;
  }

  export const groupBy = (list, fn) => {
    const groups = {};
    list.forEach(function (o) {
        const group = JSON.stringify(fn(o));
        groups[group] = groups[group] || [];
        groups[group].push(o);
    });
    return groups;
}

  export function getMaxlength(arr){
    let max = []
    arr.forEach(element => {
        max.push(TextSize.measuredTextSize(element,10).width)
    });
    return Math.max.apply(null,max)
}

export interface Node {
    name: string;
    length: number;
    children: Node[];
}

export function parse(s: any): Node {
    const ancestors = [];
    let tree = {} as Node;
    let subtree: Node;
    s = s.data
    const tokens = s.split(/\s*(;|\(|\)|,|:)\s*/);

    for (let i = 0; i < tokens.length; i++) {
        const token = tokens[i];
        switch (token) {
            case "(": // new branchset
                subtree = {} as Node;
                tree.children = [subtree];
                ancestors.push(tree);
                tree = subtree;
                break;
            case ",": // another branch
                subtree = {} as Node;
                ancestors[ancestors.length - 1].children.push(subtree);
                tree = subtree;
                break;
            case ")": // optional name next
                tree = ancestors.pop();
                break;
            case ":": // optional length next
                break;
            default:
                let x = tokens[i - 1];
                if (x === ")" || x === "(" || x === ",") {
                    tree.name = token;
                } else if (x === ":") {
                    tree.length = parseFloat(token);
                }
        }
    }
    return tree;
}


export function calcuSymbollen(symbollist){
    let lenlist = []
    symbollist.forEach(element => {
        lenlist.push(TextSize.measuredTextSize(element.replace("gene_symbol:",""),10).width)
    });
    return Math.max(...lenlist)
}