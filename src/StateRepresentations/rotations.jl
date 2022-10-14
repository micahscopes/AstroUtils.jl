function rot1Vec(vecIn, a)
    c = cos(a);
    s = sin(a);
    return SVector(vecIn[1], c*vecIn[2] + s*vecIn[3], c*vecIn[3] - s*vecIn[2])
end

function rot2Vec(vecIn, a)
    c = cos(a);
    s = sin(a);
    return SVector(c*vecIn[1] - s*vecIn[3], vecIn[2], c*vecIn[3] + s*vecIn[1]) 
end

function rot3Vec(vecIn, a)
    c = cos(a);
    s = sin(a);
    return SVector(c*vecIn[1] + s*vecIn[2], c*vecIn[2] - s*vecIn[1], vecIn[3])
end