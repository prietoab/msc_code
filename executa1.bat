@ECHO OFF
cls
SET QtdeGenes=50
SET QtdeMaxGeracoes=5000
SET QtdeIndividuos=100
SET TamanhoTorneio=2
SET ProbCrossOver=0.8
SET QtdePontosCorte=2
SET parmProbMutacao=0.3
SET IntensidadeMutacao=0.1
SET parmLambda=0.01
SET RhoMinimo=0.2

Serial_novo %QtdeGenes% %QtdeMaxGeracoes% %QtdeIndividuos% %TamanhoTorneio% %ProbCrossOver% %QtdePontosCorte% %parmProbMutacao% %IntensidadeMutacao% %parmLambda% %RhoMinimo% > Serial_Novo_parm.txt

ECHO QtdeGenes=%QtdeGenes%
ECHO QtdeMaxGeracoes=%QtdeMaxGeracoes%
ECHO QtdeIndividuos=%QtdeIndividuos%
ECHO TamanhoTorneio=%TamanhoTorneio%
ECHO ProbCrossOver=%ProbCrossOver%
ECHO QtdePontosCorte=%QtdePontosCorte%
ECHO parmProbMutacao=%parmProbMutacao%
ECHO IntensidadeMutacao=%IntensidadeMutacao%
ECHO parmLambda=%parmLambda%
ECHO RhoMinimo=%RhoMinimo%

pause