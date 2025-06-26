import './styles/App.css';
import NavBar from './components/NavBar'
import Main from './Main'
import Footer from './components/Footer'
import {useEffect, useState} from "react";

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
                setIsLoading(false);
                console.log('successfully fetched tools list')
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
            <Main toolsList={toolsList} setToolsList={setToolsList}/>
            <Footer />
        </div>
    );
}

export default App;
