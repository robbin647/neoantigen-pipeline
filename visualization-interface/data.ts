import Oviz from "crux";
import * as text_size from "crux/dist/utils/text-size";
import * as d3 from "d3";
import * as ds from "d3-sankey";
import * as dsv from "d3-dsv"
import * as TextSize from "crux/dist/utils/text-size";

export let plottree = {
    treeData:null,
}

export function plotData(_data){
    plottree.treeData = _data
    this.data.treeData = plottree.treeData

}

export function seqinfoData(_data){
}

function project(x, y) {
    var angle = (x - 90) / 180 * Math.PI, radius = y;
    return [radius * Math.cos(angle), radius * Math.sin(angle)];
}
  


export  function generatePath(node1, node2,paddingLeft){
    const path = d3.path();
    path.moveTo(node1.y + paddingLeft, node1.x);
    path.bezierCurveTo(
        (node1.y + node2.y)/2 + paddingLeft, node1.x, 
        (node1.y + node2.y)/2 + paddingLeft, node2.x, 
        node2.y + paddingLeft, node2.x
    );             
    return path.toString();
}

export function svgPathCurv(a,b,curv) {
    curv = curv ? curv : 0;
    var x1, x2, y1, y2, s, k2, controX, controY, q, l, path = '';
    s = 'M' + a.x + ',' + a.y + ' ';
    x1 = a.x; x2 = b.x
    y1 = a.y; y2 = b.y
    k2 = -(x2 - x1) / (y2 - y1);
    if(k2 < 2 && k2 > -2) {
        controX = (x2 + x1) / 2 + curv * 30;
        controX = controX < 0 ? -controX : controX;
        controY = k2 * (controX - (x1 + x2) / 2) + (y1 + y2) / 2;
        controY = controY < 0 ? -controY : controY;
    } else {
        controY = (y2 + y1) / 2 + curv * 30;
        controY = controY < 0 ? -controY : controY;
        controX = (controY - (y1 + y2) / 2) / k2 + (x1 + x2) / 2;
        controX = controX < 0 ? -controX : controX;
    }
    q = 'Q' + controX + ',' + controY + ' ';
    l = x2 + ',' + y2 + ' ';
    path = s + q + l;
    return path;
}


export function createCPath(x1, y1, x2, y2,num) {
    var path = "M" + x1 + " " + y1 + " ";
    let path2 = "M" + x1 + " " + y1 + " ";
    var cx1 = x1;
    var cy1 = (y1 + y2) * num
    var cx2 = x2;
    var cy2 = (y1 + y2) * num
    var c = "C" + cx1 + " " + cy1 + "," + cx2 + " " + cy2 + ",";
    let c2 = "C" + cx1 + " " + (y1+y2)*0.8 + ","  + cx2 + " " + (y1+y2)*0.5 + ",";
    path = path + c + x2 + " " + y2;
    path2 = path2 + c2 + x2 + " " + y2;
    return path2;
}

