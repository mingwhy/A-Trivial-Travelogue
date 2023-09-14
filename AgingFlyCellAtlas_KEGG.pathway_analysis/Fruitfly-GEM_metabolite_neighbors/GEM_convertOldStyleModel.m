
%load model
load('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat');
model=fruitflyGEM;
model.grRules{10}
model_old=convertOldStyleModel(model) %https://opencobra.github.io/cobratoolbox/stable/modules/base/utilities/index.html
model_old.rules(1)



rules=model_old.rules;
save rules.mat rules


save(model_old.rules,'rules.mat')
%%%%%%%%%%%%%%%%%%%%%
%model = creategrRulesField(model); %https://opencobra.github.io/cobratoolbox/stable/modules/base/utilities/index.html
%model.grRules{10}

names=fieldnames(model_old);
strjoin(names,'","')
%'S","annotation","b","c","compNames","comps","csense","description","eccodes","geneShortNames","genes","grRules","id","inchis","lb","metCharges","metComps","metFormulas","metNames","mets","name","osenseStr","rules","rxnConfidenceScores","rxnGeneMat","rxnNames","rxnNotes","rxnReferences","rxns","subSystems","ub","version'
%c("S","annotation","b","c","compNames","comps","csense","description","eccodes","geneShortNames","genes","grRules","id","inchis","lb","metCharges","metComps","metFormulas","metNames","mets","name","osenseStr","rules","rxnConfidenceScores","rxnGeneMat","rxnNames",
%              "rxnNotes","rxnReferences","rxns","subSystems","ub","version");

save('Fruitfly-GEM_oldStyle.mat','model')

%http://gibbs.unal.edu.co/cobradoc/cobratoolbox/tutorials/base/IO/iframe_tutorial_IO.html
%https://raw.githubusercontent.com/opencobra/COBRA.tutorials/master/base/IO/tutorial_IO.m
model_old.subSystems = vertcat(model_old.subSystems{:});
writeCbModel(model_old) %did not work
writeCbModel(model, 'format','mat', 'Fruitfly-GEM-rule.mat', 'TestModel.mat')

