class TimeHandler( object ):
    def __init__( self, timeStep, finalTime, initialTime=0.0 ):
        self.__timeStep = timeStep
        self.__finalTime = finalTime
        self.__initialTime = initialTime
        self.__currentTime = initialTime
        if self.__timeStep > self.__finalTime:
            raise Exception("Time step size cannot be greater than final time.")

    def advanceTime( self ):
        self.__currentTime += self.__timeStep

    def getTimeStep( self ):
        return self.__timeStep

    def getCurrentTime( self ):
        return self.__currentTime

    def getFinalTime( self ):
        return self.__finalTime

    def getInitialTime( self ):
        return self.__initialTime

    def isFinalTimeReached( self ):
        if self.__currentTime >= self.__finalTime:
            return False
        else:
            return True
