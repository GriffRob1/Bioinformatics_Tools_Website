import './styles/App.css';
import NavBar from './components/NavBar'
import Main from './Main'
import Footer from './components/Footer'
import {createContext, useEffect, useState} from "react";

export const AppContext = createContext(null)

function App() {

    const [toolsList, setToolsList] = useState('Loading...');
    const [isLoading, setIsLoading] = useState(true);

    //runs once when website first loads
    useEffect(() => {
        //initializes favorites for first-time users
        if (localStorage.getItem('favorites') === null) {
            localStorage.setItem('favorites', '[]')
        }

        //fetches tools list from server
        fetch('http://localhost:5000/tools-list')
            .then((response) => response.json())
            .then((data) => {
                setToolsList(data);
                setIsLoading(false)
            })
            .catch((err) => {
                setToolsList("Failed to connect");
                setIsLoading(false)
            });
    }, []);

    //waits for tools list to be fetched before loading page
    if (isLoading) return;

    return (
        <div className="app">
            <NavBar />
            <AppContext value={toolsList}>
                <Main/>
            </AppContext>
            <Footer />
        </div>
    );
}

export default App;
